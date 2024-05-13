Electric_consumption<-readRDS("data-raw/Electric_consumption_data.RDA")
Electric_consumption<-Electric_consumption[,-2]
colnames(Electric_consumption)<-c("Intensity","EE_Kitchen","EE_Laundry","EE_WH_AC")

Skin_segmentation<-readRDS("data-raw/Skin_segmentation_data.RDA")
colnames(Skin_segmentation)<-c("Skin_presence","Red","Green","Blue")

Bike_sharing<-readRDS("data-raw/Bike_sharing_data.RDA")
colnames(Bike_sharing)<-c("Rented_Bikes","Temperature","Humidity","Windspeed")

usethis::use_data(Electric_consumption,overwrite = TRUE,compress = "xz")
usethis::use_data(Skin_segmentation,overwrite = TRUE,compress = "xz")
usethis::use_data(Bike_sharing,overwrite = TRUE,compress = "xz")
