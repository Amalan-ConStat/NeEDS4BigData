Electric_consumption<-readRDS("data-raw/Electric_consumption_data.RDA")
Electric_consumption<-Electric_consumption[complete.cases(Electric_consumption),]
colnames(Electric_consumption)<-c("Intensity","Voltage","EE_Kitchen","EE_Laundry","EE_WH_AC")
Electric_consumption$Intensity<-log(Electric_consumption$Intensity)

Skin_segmentation<-readRDS("data-raw/Skin_segmentation_data.RDA")
colnames(Skin_segmentation)<-c("Skin_presence","Red","Green","Blue")

Bike_sharing<-readRDS("data-raw/Bike_sharing_data.RDA")
colnames(Bike_sharing)<-c("Rented_Bikes","Temperature","Humidity","Windspeed")

One_Million_Songs<-readRDS("data-raw/One_Million_Songs_data.RDA")
colnames(One_Million_Songs)<-c("Counts","Duration","Loudness",
                               "Tempo","Artist_Hotness","Song_Hotness")

usethis::use_data(Electric_consumption,overwrite = TRUE,compress = "xz")
usethis::use_data(Skin_segmentation,overwrite = TRUE,compress = "xz")
usethis::use_data(Bike_sharing,overwrite = TRUE,compress = "xz")
usethis::use_data(One_Million_Songs,overwrite = TRUE,compress = "xz")
