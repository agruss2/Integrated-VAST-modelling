###############################################################################################################################
##
##  (04) This script runs a model fitted to COD BT data for HAK4
##  
###############################################################################################################################

######## Set the output directory
if ( !dir.exists( paste0( DIR$Output, "\\HAK4_COD_BT_data" ) ) ) {
	dir.create( paste0( DIR$Output, "\\HAK4_COD_BT_data" ) )
} 
OutputDir <- paste0( DIR$Output, "\\HAK4_COD_BT_data" )
setwd( OutputDir )

######## Load required R packages
library( VAST )
library( DHARMa )

######## Load HAK data and keep data only for HAK4
load( make.filename( "HAK_Sp_data.RData", DIR$Input ) )
Sp_data <- Sp_data[Sp_data$Dataset == "COD_BT",]
dim( Sp_data )
names( Sp_data )
####  [1] "Lon"                "Lat"                "Year"              
####  [4] "AreaSwept_km2"      "Value"              "Dataset"           
####  [7] "Datatype"           "Vessel_ID"          "Monitoring_program"
#### [10] "QMA"      

nrow( Sp_data ) #### 150 928 data points
Sp_data <- Sp_data[Sp_data$Year > 1990,]
Sp_data <- Sp_data[Sp_data$QMA == "HAK4",]
nrow( Sp_data ) #### 53 783 data points
Sp_data$Monitoring_program <- droplevels( Sp_data$Monitoring_program )
Sp_data$Vessel_ID <- droplevels( Sp_data$Vessel_ID )
table( Sp_data$Vessel_ID )
Sp_data$Year <- as.factor( Sp_data$Year )
Sp_data$Year <- droplevels( Sp_data$Year )
Sp_data$Year <- as.numeric( as.character( Sp_data$Year ) )
table( Sp_data$Year )
#### 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 
#### 1132  607 1113 1533 1356  780 1132 1650 1776 1898 1635 2162 1357 1348 1680 1658 
#### 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 
#### 2006 2315 2394 1619 1612 1841 1777 1263 2317 1576 1925 1428 2124 2123 2733 1913

years <- sort( unique( Sp_data$Year ) )
Year_i <- as.numeric( as.character( Sp_data$Year ) )
Lon_i <- Sp_data$Lon
Lat_i <- Sp_data$Lat
table( Sp_data$Dataset )
table( Sp_data$Datatype )
tapply( Sp_data[,'Value'], INDEX = list( Sp_data[,'Year'], Sp_data[,'Dataset'] ), FUN = mean )
table( ifelse( Sp_data$Value > 0, 1, 0 ), Sp_data$Year )
####     1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005
####   0  496  239  529  967 1017  275  745  502  690  979  806  803  901  811  956
####   1  636  368  584  566  339  505  387 1148 1086  919  829 1359  456  537  724
####    
####     2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020
####   0  939 1202 1513 1669 1137  863 1001  725  519 1456  947 1077  922 1466 1340
####   1  719  804  802  725  482  749  840 1052  744  861  629  848  506  658  783
####    
####     2021 2022
####   0 1698 1309
####   1 1035  604

######## Load the prediction grid for HAK4
load( make.filename( "Prediction_grid_HAK4_10kmx10km.RData", DIR$Input ) )
dim( Prediction_grid )
names( Prediction_grid )
input_grid <- as.data.frame( cbind( Prediction_grid$Lon, Prediction_grid$Lat, Prediction_grid$Area_km2 ) )
names( input_grid ) <- c( "Lon", "Lat", "Area_km2" )
summary( input_grid$Lon )
summary( input_grid$Lat )

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
strata.limits <- data.frame( STRATA = "All_areas" )

######## Define the extrapolation grid
Extrapolation_List = make_extrapolation_info( Region = Region, strata.limits = strata.limits, 
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
Q_ik = NULL
v_i <- as.numeric( as.factor( Sp_data$Vessel_ID ) )
TmbData = make_data( "Version" = Version, "FieldConfig" = FieldConfig, "OverdispersionConfig" = OverdispersionConfig,
	"RhoConfig" = RhoConfig, "ObsModel" = ObsModel, "c_i" = rep( 0, nrow( Sp_data ) ), 
	"b_i" = as.numeric( as.character( Sp_data[, "Value"] ) ), 
	"a_i" = rep( 1, nrow( Sp_data ) ), 
	"v_i" = v_i - 1,
	"s_i" = Sp_data[, "knot_i"] - 1, 
	"t_i" = as.numeric( as.character( Sp_data[, "Year"] ) ),
	"a_xl" = Spatial_List[["a_gl"]], "Q_ik" = Q_ik, "MeshList" = Spatial_List$MeshList,
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

