library(sf)
library(stats19)
library(dplyr)
library(tidyr)
library(ggplot2)
library(osmextract)
library(osmdata)
library(tmap)
library(spdep)
library(sfdep)
library(pct)
library(readr)

#importing the datasets
?get_stats19

ac_2018 = get_stats19(year=2018,type = "ac")
ac_2019 = get_stats19(year=2019,type = "ac")
ac_2020 = get_stats19(year=2020,type = "ac")
ac_2021 = get_stats19(year=2021,type = "ac")
ac_2022 = get_stats19(year=2022,type = "ac")


cas_2018 = get_stats19(year=2018,type = "cas")
cas_2019 = get_stats19(year=2019,type = "cas")
cas_2020 = get_stats19(year=2020,type = "cas")
cas_2021 = get_stats19(year=2021,type = "cas")
cas_2022 = get_stats19(year=2022,type = "cas")

write_csv(cas_2018,"data/cas_2018.csv")
write_csv(cas_2019,"data/cas_2019.csv")
write_csv(cas_2020,"data/cas_2020.csv")
write_csv(cas_2021,"data/cas_2021.csv")
write_csv(cas_2022,"data/cas_2022.csv")

write_csv(ac_2018,"data/ac_2018.csv")
write_csv(ac_2019,"data/ac_2019.csv")
write_csv(ac_2020,"data/ac_2020.csv")
write_csv(ac_2021,"data/ac_2021.csv")
write_csv(ac_2022,"data/ac_2022.csv")
#pre-processing

colnames(ac_2019) == colnames(ac_2022)
colnames(cas_2019) == colnames(cas_2022)

unique(cas_2021$casualty_class)
colnames(cas_2021)
cas_2019$lsoa_of_casualty = as.character(cas_2019$lsoa_of_casualty)

ac_comb = bind_rows(ac_2018,ac_2019,ac_2020,ac_2021,ac_2022)
cas_comb = bind_rows(cas_2018,cas_2019,cas_2020,cas_2021,cas_2022)

# null and inf values

sum(is.na(ac_comb))
sum(is.na(cas_comb))

sum(is.infinite(ac_comb$accident_index))
sum(is.infinite(cas_comb$accident_index))


# hotspot analysis by lsoa

colnames(ac_comb)

lsoa_crash = ac_comb %>% 
  filter(local_authority_district == "Leeds") %>% 
  select("lsoa_of_accident_location",
         "number_of_casualties") %>% 
  group_by(lsoa_of_accident_location) %>% 
  summarize(casualties = sum(number_of_casualties, na.rm=TRUE),
            accidents = n())

#getting external data for lsoa boundaries

unzip("Lower_layer_Super_Output_Areas_2021_EW_BFC_V8_8154990398368723939.zip")
lsoa = read_sf("LSOA_2021_EW_BFC_V8.shp")
unique(lsoa$LSOA21NM)

tinytex::install_tinytex()
install.packages("tinytex")

leeds_lsoa = lsoa[grep("Leeds", lsoa$LSOA21NM), ]
nrow(leeds_rows)
st_write(leeds_lsoa,"leeds_lsoa.shp")

leeds_lsoa = read_sf("leeds_lsoa.shp")
accidents = leeds_lsoa %>% 
  left_join(lsoa_crash, by = c("LSOA21CD"="lsoa_of_accident_location"))

head(accidents)

DataExplorer::plot_missing(accidents)

accidents = accidents %>%
  mutate(accidents = if_else(is.na(accidents), 0, accidents),
         casualties = if_else(is.na(casualties), 0, casualties))

tm_shape(accidents) +
  tm_fill(col = "casualties", style = "quantile", palette = "Blues") +
  tm_borders() +
  tm_layout(legend.position = c("left", "bottom"))

tm_shape(accidents) +
  tm_fill(col = "accidents", style = "quantile", palette = "Blues") +
  tm_borders() +
  tm_layout(legend.position = c("left", "bottom"))


#hotspot analysis using Getis ord Gi

ac_nb = poly2nb(accidents, queen = TRUE)
View(ac_nb)

#check for empty neighbours
empty_nb <- which(card(ac_nb) == 0)
empty_nb  

# Binary weighting assigns a weight of 1 to all neighboring features and a weight of 0 to all other features
ac_w_binary = nb2listw(ac_nb, style="B")

# Calculate spatial lag of casualties
ac_lag = lag.listw(ac_w_binary, accidents$casualties)  

globalG.test(accidents$casualties, ac_w_binary)

# Identify neighbors, create weights, calculate spatial lag
ac_nbs = accidents |> 
  mutate(
    nb = st_contiguity(geometry),        # neighbors share border/vertex
    wt = st_weights(nb),                 # row-standardized weights
    ac_lag = st_lag(casualties, nb, wt)    # calculate spatial lag of TreEqty
  ) 


# Calculate the Gi using local_g_perm
ac_hot_spots = ac_nbs |> 
  mutate(
    Gi = local_g_perm(casualties, nb, wt, nsim = 999)
    # nsim = number of Monte Carlo simulations (999 is default)
  ) |> 
  # The new 'Gi' column itself contains a dataframe 
  # We can't work with that, so we need to 'unnest' it
  unnest(Gi) 

ac_hot_spots |> 
  ggplot((aes(fill = gi))) +
  geom_sf(color = "black", lwd = 0.15) +
  scale_fill_gradient2() # makes the value 0 (random) be the middle

# Create a new data frame called hpa
hpa = ac_hot_spots |> 
  # with the columns 'gi' and 'p_folded_sim"
  # 'p_folded_sim' is the p-value of a folded permutation test
  select(LSOA21CD, gi, p_folded_sim) |> 
  mutate(
    # Add a new column called "classification"
    classification = case_when(
      # Classify based on the following criteria:
      gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
      gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
      gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
      gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
      gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
      gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
      TRUE ~ "Insignificant"
    ),
    # Convert 'classification' into a factor for easier plotting
    classification = factor(
      classification,
      levels = c("Very hot", "Hot", "Somewhat hot",
                 "Insignificant",
                 "Somewhat cold", "Cold", "Very cold")
    )
  ) 

# Visualize the results with ggplot2
ggplot(data= hpa, aes(fill = classification)) +
geom_sf(color = "black", lwd = 0.1) +
scale_fill_brewer(type = "div", palette = 5) +
theme_void() +
labs(
  fill = "Hot Spot Classification",
  title = "Hotspot Analysis by Casualties"
)

tm_shape(hpa) +
  tm_fill(col = "classification", style = "quantile") +
  tm_borders() +
  tm_layout(legend.position = c("left", "bottom"))



#road data analysis

sf_osm_leeds = 
  opq("leeds uk") %>% 
  add_osm_feature(key = "highway", value = c("motorway",
                                             "trunk",
                                             "primary",
                                             "secondary",
                                             "tertiary",
                                             "residential",
                                             "unclassified",
                                             "road")) %>% 
  osmdata_sf()


road_lds = sf_osm_leeds$osm_lines %>% select("osm_id", "name", "highway", "geometry")

sf_osm_leeds_multiline <-
  road_lds %>%
  filter(!is.na(name)) %>% 
  group_by(name) %>% 
  summarise(count = n()) %>% 
  mutate(name = as.character(name))

crash = ac_comb %>% 
  filter(local_authority_district == "Leeds") %>% 
  select("lsoa_of_accident_location",
         "number_of_casualties",
         "latitude",
         "longitude",
         "location_easting_osgr",
         "location_northing_osgr",
         "speed_limit")
  
crash_sf = st_as_sf(crash,coords = c("location_easting_osgr","location_northing_osgr"),crs="EPSG:27700")
roads_lds = st_transform(sf_osm_leeds_multiline,crs="EPSG:27700")

road_crashes = roads_lds %>% 
  st_join(crash_sf, join = st_is_within_distance, dist = 100) %>%
  select("name","number_of_casualties","speed_limit") %>%  
  group_by(name) %>% 
  summarize(n=n(), number_of_casualties = sum(number_of_casualties),speed_limit=mean(speed_limit))%>% 
  filter(!is.na(name)) %>% 
  arrange(desc(n))

plot(road_crashes[2])

a = road_crashes %>% st_drop_geometry()
plot(as.list(a$n),as.list(a$number_of_casualties))

# traffic lights

sf_tr = 
  opq("leeds uk") %>% 
  add_osm_feature(key = "highway", value = c("traffic_signals")) %>% 
  osmdata_sf()

trfc_lts = st_transform(sf_tr$osm_points,crs="EPSG:27700")

plot(trfc_lts$geometry)

traffic_lights = trfc_lts %>% 
  st_join(roads_lds, join = st_is_within_distance, dist = 50) %>%
  select("name.y") %>%  
  group_by(name.y) %>% 
  summarize(n_tl=n())%>% 
  rename(name = name.y) %>% 
  filter(!is.na(name)) %>% 
  arrange(desc(n_tl))

traffic_lights

#traffic calming measures

sf_tr_cal = 
  opq("leeds uk") %>% 
  add_osm_feature(key = "traffic_calming", value = c("bump",
                                                     "mini_bumps",
                                                     "hump",
                                                     "table",
                                                     "cushion")) %>% 
  osmdata_sf()

tr_cal = st_transform(sf_tr_cal$osm_points, crs="EPSG:27700")

traffic_calming = tr_cal %>% 
  st_join(roads_lds, join = st_is_within_distance, dist = 50) %>%
  select("name") %>%  
  group_by(name) %>% 
  summarize(n_tc=n()) %>%
  filter(!is.na(name)) %>% 
  arrange(desc(n_tc))

traffic_calming


#nightclub

sf_nc = 
  opq("leeds uk") %>% 
  add_osm_feature(key = "amenity", value = c("bar",
                                             "pub",
                                             "nightclub")) %>% 
  osmdata_sf()

nc = st_transform(sf_nc$osm_points, crs="EPSG:27700")

night_clubs = nc %>% 
  st_join(roads_lds, join = st_is_within_distance, dist = 500) %>%
  select("name.y") %>%  
  group_by(name.y) %>%
  rename(name=name.y) %>% 
  summarize(n_nc=n()) %>%
  filter(!is.na(name)) %>% 
  arrange(desc(n_nc))

night_clubs




# regression

attr_final = st_drop_geometry(road_crashes) %>% 
  left_join(st_drop_geometry(traffic_lights),by="name") %>% 
  left_join(st_drop_geometry(traffic_calming, by="name")) %>% 
  left_join(st_drop_geometry(night_clubs, by="name"))

attr_final[is.na(attr_final)] = 0
attr_final

lr_model = lm(n~n_tl+n_tc+n_nc+speed_limit, data=attr_final)
summary(lr_model)

colnames(ac_comb)



# pct

zones=get_pct(region = "west-yorkshire",geography="lsoa",layer="l")
unique(zones$geo_name1)
leeds_dl = zones[grep("Leeds", zones$geo_name1), ]
plot(st_transform(leeds_dl$geometry,crs="EPSG:27700"),add=TRUE)
leeds_dl %>% arrange(desc(all)) %>% select("all","bicycle","foot","car_driver")
