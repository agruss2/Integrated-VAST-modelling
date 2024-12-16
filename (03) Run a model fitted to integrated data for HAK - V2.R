###############################################################################################################################
##
##  (03) This script runs a model fitted to integrated data (i.e., TRAWL BT data + COD BT data) for HAK -
##  Version 2 (V2): Here the model is a spatially-varying catchability model
##  
###############################################################################################################################

######## Set the output directory
if ( !dir.exists( paste0( DIR$Output, "\\HAK_Integrated_data_V2" ) ) ) {
	dir.create( paste0( DIR$Output, "\\HAK_Integrated_data_V2" ) )
} 
OutputDir <- paste0( DIR$Output, "\\HAK_Integrated_data_V2" )
setwd( OutputDir )

######## Load required R packages
library( VAST )
library( DHARMa )

######## Load HAK data 
load( make.filename( "HAK_Sp_data.RData", DIR$Input ) )
dim( Sp_data )
names( Sp_data )
####  [1] "Lon"                "Lat"                "Year"              
####  [4] "AreaSwept_km2"      "Value"              "Dataset"           
####  [7] "Datatype"           "Vessel_ID"          "Monitoring_program"
#### [10] "QMA"      

nrow( Sp_data ) #### 156 980 data points
Sp_data_1 <- Sp_data[Sp_data$Dataset == "TRAWL_BT",]
Sp_data_1 <- Sp_data_1[Sp_data_1$Monitoring_program %in% c( "CHAT_MD", "SOUTH_MD", "SUBA_AUT", "SUBA_SUM", "WCSI_MD" ),]
Sp_data_2 <- Sp_data[Sp_data$Dataset == "COD_BT",]
Sp_data_2 <- Sp_data_2[Sp_data_2$Year > 1990,] 
Sp_data <- rbind( Sp_data_1, Sp_data_2 )
nrow( Sp_data ) #### 140 229 data points
Sp_data$Monitoring_program <- droplevels( Sp_data$Monitoring_program )
Sp_data$Vessel_ID <- droplevels( Sp_data$Vessel_ID )
table( Sp_data$Vessel_ID )
Sp_data$Year <- as.factor( Sp_data$Year )
Sp_data$Year <- droplevels( Sp_data$Year )
Sp_data$Year <- as.numeric( as.character( Sp_data$Year ) )
table( Sp_data$Year )
#### 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 
#### 2577 2969 3632 3160 2415 1872 2320 3317 3940 4193 4612 4851 3745 3392 3905 4175 
#### 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 
#### 4652 5389 4724 4094 4408 4481 5703 4383 4860 4656 6101 6539 6220 7056 7137 4751 

years <- sort( unique( Sp_data$Year ) )
Year_i <- as.numeric( as.character( Sp_data$Year ) )
Lon_i <- Sp_data$Lon
Lat_i <- Sp_data$Lat
table( Sp_data$Dataset )
table( Sp_data$Datatype )
tapply( Sp_data[,'Value'], INDEX = list( Sp_data[,'Year'], Sp_data[,'Dataset'] ), FUN = mean )
table( ifelse( Sp_data$Value > 0, 1, 0 ), Sp_data$Year )
####     1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005
####   0 1320 1439 1696 2085 1647  748 1558 1514 1714 1817 2549 2211 2231 2223 2614
####   1 1257 1530 1936 1075  768 1124  762 1803 2226 2376 2063 2640 1514 1169 1291
####    
####     2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020
####  0 2558 3016 3629 3255 2641 2758 2598 3011 2263 2989 2983 3602 3669 4396 4894
####   1 1617 1636 1760 1469 1453 1650 1883 2692 2120 1871 1673 2499 2870 1824 2162
####    
####     2021 2022
####   0 4689 3149
####   1 2448 1602

######## Load the prediction grid for HAK
load( make.filename( "Prediction_grid_HAK_10kmx10km.RData", DIR$Input ) )
dim( Prediction_grid )
names( Prediction_grid )
table( Prediction_grid$QMA )
#### HAK1 HAK4 HAK7 
#### 7664 3995 2393 

input_grid <- as.data.frame( cbind( Prediction_grid$Lon, Prediction_grid$Lat, 
	Prediction_grid$Area_km2, as.factor( Prediction_grid$QMA ) ) )
names( input_grid ) <- c( "Lon", "Lat", "Area_km2", "strata" )
summary( input_grid$Lon )
summary( input_grid$Lat )
table( input_grid$strata )

######## Define some VAST settings 
Version = get_latest_version( package = "VAST" )
Method = "Barrier"
grid_size_km = 25
n_x = 200
FieldConfig = matrix( data = c( Omega1 = "IID", Epsilon1 = "IID", Beta1 = "IID",  
	Omega2 = "IID", Epsilon2 = "IID", Beta2 = "IID" ), nrow = 3, ncol = 2 )
RhoConfig = c( Beta1 = 0, Beta2 = 0, Epsilon1 = 4, Epsilon2 = 4 ) 
OverdispersionConfig = c( Eta1 = 1, Eta2 = 1 )
ObsModel = c( 2, 1 )
fine_scale = TRUE

######## Decide which post-hoc calculations to include in VAST output
Options =  c( "SD_site_density" = 1, "SD_site_logdensity" = 0, "Calculate_Range" = 1, "Calculate_evenness" = 0, 
	"Calculate_effective_area" = 1, "Calculate_Cov_SE" = 0, "Calculate_Synchrony" = 0, 
	"Calculate_Coherence" = 0, "report_additional_variables" = TRUE, 	
	"Range_fraction" = 0.2 )

######## Determine the study region 
Region = "User"

######## Determine strata within the study region 
strata.limits = data.frame( 'STRATA' = sort( unique( input_grid$strata ) ) )

######## Define the extrapolation grid
Extrapolation_List = make_extrapolation_information( Region = Region, strata.limits = strata.limits, 
	input_grid = input_grid, max_cells = 2000 ) 

######## Generate the spatial information necessary for conducting spatio-temporal parameter estimation
Spatial_List = make_spatial_info( grid_size_km = grid_size_km, n_x = n_x, fine_scale = fine_scale, Method = Method, 
	Lon = as.numeric( as.character( Sp_data[,"Lon"] ) ),
	Lat = as.numeric( as.character( Sp_data[, "Lat"] ) ), 
	Extrapolation_List = Extrapolation_List, Save_Results = TRUE, "knot_method" = "grid" )
Sp_data = cbind( Sp_data, knot_i = Spatial_List$knot_i )
save( Sp_data, file = "Sp_data.RData" )

######## Plot data and knots 
Sp_data$Year <- as.numeric( Sp_data$Year )
plot_data( Extrapolation_List = Extrapolation_List, Spatial_List = Spatial_List, Data_Geostat = Sp_data ) 

######## Build the TMB object
Sp_data$Type <- "5" 
Sp_data$Type <- ifelse( Sp_data$Monitoring_program == "CHAT_MD", "0", Sp_data$Type ) 
Sp_data$Type <- ifelse( Sp_data$Monitoring_program == "SUBA_SUM", "1", Sp_data$Type ) 
Sp_data$Type <- ifelse( Sp_data$Monitoring_program == "SOUTH_MD", "2", Sp_data$Type ) 
Sp_data$Type <- ifelse( Sp_data$Monitoring_program == "WCSI_MD", "3", Sp_data$Type ) 
Sp_data$Type <- ifelse( Sp_data$Monitoring_program == "SUBA_AUT", "4", Sp_data$Type ) 
Sp_data$Type <- as.factor( Sp_data$Type )
table( Sp_data$Type )
Q1_formula <- ~ Type
Q1config_k <- c( 3, 3, 3, 3, 3 ) 
Q2_formula <- ~ 1
Q2config_k <- NULL
catchability_data <- Sp_data[,-8]
v_i <- as.numeric( as.factor( Sp_data$Vessel_ID ) )
TmbData = make_data( "Version" = Version, "FieldConfig" = FieldConfig, "OverdispersionConfig" = OverdispersionConfig,
	"RhoConfig" = RhoConfig, "ObsModel" = ObsModel, "c_i" = rep( 0, nrow( Sp_data ) ), 
	"b_i" = as.numeric( as.character( Sp_data[, "Value"] ) ), 
	"a_i" = rep( 1, nrow( Sp_data ) ), 
	"v_i" = v_i - 1, 
	"s_i" = Sp_data[, "knot_i"] - 1, 
	"t_i" = as.numeric( as.character( Sp_data[, "Year"] ) ),
	"a_xl" = Spatial_List[["a_gl"]], "catchability_data" = catchability_data,
  	"Q1_formula" = Q1_formula, "Q2_formula" = Q2_formula, "Q1config_k" = Q1config_k, "Q2config_k" = Q2config_k,
	"MeshList" = Spatial_List$MeshList,
	"GridList" = Spatial_List$GridList, "Method" = Spatial_List$Method, "Options" = Options, 
	"CheckForErrors" = TRUE, "spatial_list" = Spatial_List )

######## Build the VAST model 
TmbList = make_model( "build_model" = TRUE, "TmbData" = TmbData, "RunDir" = OutputDir, "Version" = Version, 
	"RhoConfig" = RhoConfig, "loc_x" = Spatial_List$loc_x, "Method" = Method, "Use_REML" = FALSE )

######## Check parameters
Obj = TmbList[["Obj"]]
Obj$fn( Obj$par )
Obj$gr( Obj$par )
any( Obj$gr( Obj$par ) == 0 )

######## Estimate fixed effects and predict random effects
Opt = TMBhelper::fit_tmb( obj = Obj, lower = TmbList[["Lower"]], upper = TmbList[["Upper"]], 
	getsd = TRUE, savedir = OutputDir, bias.correct = TRUE, newtonsteps = 3, 
	bias.correct.control = list( sd = FALSE, split = NULL, nsplit = 1, vars_to_correct = c( "Index_cyl", "Index_ctl" ) ), 
	getJointPrecision = TRUE ) 
Report = Obj$report()
Save = list( "Opt" = Opt, "Report" = Report, "ParHat" = Obj$env$parList( Opt$par ), "TmbData" = TmbData )
save( Save, file = "Save.RData" )
save( Obj, file = "Obj.RData" )

######## Print the diagnostics generated during parameter estimation, and confirm that:
######## (1) no parameter is hitting an upper or lower bound and (2) the final gradient for each fixed-effect 
######## is close to zero (less than 0.0001). Also check model convergence via the Hessian (should be TRUE)
pander::pandoc.table( Opt$diagnostics[,c( 'Param', 'Lower', 'MLE', 'Upper', 'final_gradient' )] ) 
all( abs( Opt$diagnostics[,'final_gradient'] ) <1e-4 ) #### TRUE
all( eigen( Opt$SD$cov.fixed )$values >0 ) #### TRUE

