library(tidyverse)
library(sf)
library(did)
library(rnaturalearth)
library(rnaturalearthdata)
library(geodata) 
library(exactextractr) 


first_year = 1980
last_year = 2019
treatment_radius = 25000 


##### Download geodata

population_data <- geodata::population(year = 2000, res = 0.5, path = "pop_data") # Download population raster from geodata 

canada_outline <- ne_countries(scale = "large", returnclass = "sf") |> # Canada outline from natural earth
  filter(name == "Canada")

##### Prepare conflict data

conflicts_raw <- read.csv("data/gdelt_canada_full.csv") # GDELT V1 full data for Canada

conflicts <- conflicts_raw |> 
  filter(ActionGeo_FullName != "Canada") # Filter out data which GDELT dosen't know the location of

protests <- conflicts |> 
  filter(EventRootCode == 14) # Filter only for protest events


protests <- protests |> 
  st_as_sf(coords = c("ActionGeo_Long", "ActionGeo_Lat"), crs = 4326) # Turn into spatial data frame, for analysis

protests <- protests |> 
  st_filter(canada_outline) # Filter out events which are geolocated outside of Canada 

protests <- protests |> 
  mutate(
    date = as.Date(as.character(SQLDATE), format = "%Y%m%d"), # convert SQLDATE first to string then to date with specific format
  ) |> 
  mutate(year = year(date))  # extract year into own column with year function


protests <- protests |> 
  filter(year != 1920) |> # some protests are falsely labeled as year 1920 
  filter(year >= first_year) |> # select for year range
  filter(year <= last_year)

protests |> 
  st_drop_geometry() |> 
  group_by(year) |> 
  summarise(count = n()) |> 
  ungroup() |> 
  ggplot() +
  geom_col(aes(x = year, y = count)) +
  labs(
    title = "Figure 1. Protest Events per Year"
  )


protests <- protests |> 
  select(-SQLDATE, -ActionGeo_FullName, -EventRootCode, -date) # Remove columns irrelevant for did: keep year, geodata

protest_places <- protests |> # Prepare places of protests for later mapping
  group_by(geometry) |> 
  summarise(count = n()) |> 
  ungroup() |> 
  arrange(desc(count)) 


##### Prepare Dam data 

dams <- read_sf("data/GDW_v1_0_shp/GDW_barriers_v1_0.shp") |> 
  filter(COUNTRY == "Canada") |>  # filter for Canada
  mutate(year = YEAR_DAM)

dams <- dams |> 
  filter(year >= first_year) |> 
  filter(year <= last_year)


# Plot dams and protests
ggplot() +
  geom_sf(data = canada_outline) +
  geom_sf(data = protest_places ,size = 0.6, aes(color = "Protest")) + 
  geom_sf(data = dams, size = 0.6, aes(color = "Dam")) +
  scale_color_manual(
    name = "",
    values = c("Protest" = "red", "Dam" = "blue")
  ) +
  labs(
    title = "Figure 2. Protests and Dams in Canada"
  )



##### Prepare data for did

crs <- 3978 # I use this north america specific coordinate reference system because lat/long is not suitable for calculating distances -> projects degrees into meters

 
### Select relevant columns for did (in this case only year and location). Transform to correct crs. Attach unique ids
dams <- dams |>  
  select(year) |> 
  st_transform(crs = crs) |> 
  mutate(dam_id = row_number())

protests <- protests |> 
  select(year) |>
  st_transform(crs = crs) |> 
  mutate(protest_id = row_number())


canada_outline <- canada_outline |> 
  st_transform(crs = crs)

dam_buffers <- st_buffer(dams, dist = treatment_radius) |>  # Create circles around dams to act as treatment zones
  mutate(unit_id = paste0("dam_", dam_id), # give unique id to the zone, called unit_id
         first_treat_year = (year - 2), # first_treat_year is the year the dam was built - 2 to account for protests before
         is_treated = 1) |> # indicate as treatment zone
  select(unit_id, first_treat_year, is_treated) # select for did relevant columns

circle_union <- st_union(dam_buffers) # needed to know where not to make grid, st_filter only works with one object

grid <- st_make_grid(canada_outline, # make grid
                     cellsize = treatment_radius * 2, # grids should be roughly same size as treatment zones
                     square = TRUE) |> 
  st_intersection(canada_outline) # filter grid to be only within Canada

grid <- st_sf(cell_id = 1:length(grid), geometry = grid) # make data frame out of grid and add cell_id

control_cells <- grid |>
  st_filter(circle_union, .predicate = st_disjoint) |> # filter out all cells which intersect a dam (treatment zone)
  mutate(unit_id = paste0("control_", cell_id), # give unique id just like dam circles
         first_treat_year = 0, # never treated = 0
         is_treated = 0) |> # mark as control group 
  select(unit_id, first_treat_year, is_treated) # select for did relevant columns

all_cells <- bind_rows(dam_buffers, control_cells) # make one data frame for all cells treatment and control, works because they have the same columns

#### Calculate population for all cells

all_cells <- all_cells |> 
  st_transform(crs = 4326) # convert back to wgs84 because population data is in this format

all_cells$population <- exact_extract(population_data, all_cells, 'sum') # for each cell sum up all population within

# Population of all cells 
ggplot() +
  geom_sf(data = canada_outline) +
  geom_sf(data = all_cells, aes(fill = population)) +
  scale_fill_continuous(name = "Population", label = scales::comma) +
  labs(
    title = "All Units with their Population",
    caption = "Grid cells are control units and circles are treatment units."
    
  )

ggplot() +
  geom_sf(data = canada_outline) +
  geom_sf(data = (all_cells |> filter(population >= 100)), aes(fill = population)) +
  scale_fill_continuous(name = "Population", label = scales::comma) +
  labs(
    title = "Figure 3. All Units with a Population Higher than 100",
    caption = "Grid cells are control units and circles are treatment units."
  )


all_cells <- all_cells |> 
  st_transform(crs = crs) |> # convert back to Canada CRS
  filter(population >= 100) # filter out all cells with less than 100 people

all_cells |> 
  filter(is_treated == 1) |> 
  count() # count how many dams remain

all_cells |> 
  filter(is_treated == 0) |> 
  count() # count how many control groups remain


protests_matched <- st_join(
  protests,
  all_cells,
  join = st_intersects # join the cell in which the protests is located to the row of the protest 
  ) |> 
  filter(!is.na(unit_id)) # filter out all protests which are not within cells (population < 100 / space between circle and cell)


protests_matched |> 
  count() # how many protests get used in the final result


protests_counts <- protests_matched |>
  st_drop_geometry() |>             
  group_by(unit_id, year, population) |> # prepare for did and group by unit and year (and population but it dosent form more groups)
  summarise(protests_count = n(), .groups = "drop") |> # count how many protests per year per cell
  mutate(cpc = (protests_count / population) * 100000) # calculate the protest count per 100k population
# years without protests don't show up and will be added in the next step

all_years <- min(protests$year, na.rm = TRUE):max(protests$year, na.rm = TRUE) # all years in the protest data

df_panel <- expand_grid(unit_id = all_cells$unit_id, year = all_years) |> # build data frame with all combinations of unit_id and year, because summarize removed 0 protest years I have to use expand_grid function
  left_join(all_cells |> select(unit_id, first_treat_year, is_treated, population) |> st_drop_geometry(), by = "unit_id") |> # join back the information about the units which was lost due to expandgrid
  left_join(protests_counts |> select(unit_id, year, cpc), by = c("unit_id", "year")) |> # join back cpc
  mutate(cpc = coalesce(cpc, 0)) # replace NA's in cpc with 0s

df_panel <- df_panel |>
  mutate(unit_id_num = as.integer(factor(unit_id))) # make unit id numeric with factor function which gives each unit id a unique number

df_panel <- df_panel |> 
  mutate(log_cpc = log1p(cpc)) # make cpc logarithmic with ln(cpc + 1), + 1 is necessary otherwise 0 turns to -inf

# call att_gt from did package 
result <- att_gt(
  yname = "log_cpc", # outcome -> logged cpc
  tname = "year", # time variable
  idname = "unit_id_num", # cell id 
  gname = "first_treat_year", # treatment group -> dam build year
  data = df_panel,
  control_group = "nevertreated", # uses never treated groups as control groups
  clustervars = "unit_id_num", #how data should be grouped -> after dam units
  est_method = "dr",
  xformla = ~ log(population) # population is my only covariate
)

summary(result)

agg_simple <- aggte(result, type = "simple") 
summary(agg_simple)
capture.output(summary(agg_simple), file = "simple_summary.txt")


agg_dynamic <- aggte(result, type = "dynamic")
summary(agg_dynamic)
capture.output(summary(agg_dynamic), file = "dynamic_summary.txt")


ggdid(agg_dynamic) +
  ggtitle("Figure 4. Average Effect by Length of Exposure") +
  scale_x_continuous(breaks = seq(-32, 38, by = 2))


#### Appendix

dams_appendix <- read_sf("data/GDW_v1_0_shp/GDW_barriers_v1_0.shp") |> 
  filter(COUNTRY == "Canada") |>  
  mutate(year = YEAR_DAM) |> 
  filter(year >= first_year) |> 
  filter(year <= last_year) |> 
  mutate(dam_id = row_number()) |> 
  relocate(dam_id)


dams_all_table <- dams_appendix |> 
  st_drop_geometry() |> 
  rename(main_use = MAIN_USE) |> 
  arrange(year) |> 
  select(dam_id, year, main_use)


dams_treatment_table <- dams_appendix |> 
  st_transform(crs = crs) |> 
  st_filter(all_cells) |> 
  st_drop_geometry() |> 
  rename(main_use = MAIN_USE) |> 
  arrange(year) |> 
  select(dam_id, year, main_use)

dams_all_table |> write.csv("../Abgabe/plots_tables/table_3.csv")
dams_treatment_table |> write.csv("../Abgabe/plots_tables/table_4.csv")




