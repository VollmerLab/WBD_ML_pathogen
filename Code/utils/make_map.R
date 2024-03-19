library(sf)
library(tidyverse)
library(magrittr)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggstar)
library(ggmagnify)

#### Functions ####
convert_to_decimal <- function(x){
  out <- str_split(x, 'W|N') %>%
    unlist(recursive = FALSE) %>%
    as.numeric() %>%
    magrittr::divide_by(c(1, 60)) %>%
    sum
  
  if_else(str_detect(x, 'W'), -1 * out, out)
}

#### Data ####
#https://stridata-si.opendata.arcgis.com/content/4239c3bf7d534e758f66e155d0e64c5c/about
bocas_shape <- list.files(path = '../../Data/bocas_shape/commondata/', 
           pattern = '.shp$', recursive = TRUE) %>%
  tibble(basename = .) %>%
  mutate(file = str_c('../../Data/bocas_shape/commondata/', basename),
         folder = str_extract(basename, '^.*/') %>% str_remove('/'),
         basename = str_remove(basename, '\\.shp') %>% str_remove('^.*/')) %>%
  rowwise(folder, basename) %>%
  reframe(st_read(file, quiet = TRUE)) %>%
  st_as_sf()


world <- ne_countries(scale = "large", returnclass = "sf")

site_list <- read_csv('C:/Users/jdsel/Dropbox/1_AC_field_abundprev/Acer5sites_HS_Added.csv', 
         show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  select(site) %>%
  distinct

panama_locations <- tribble(
  ~'reef', ~'lon_garmin', ~'lat_garmin', 
  'Tetas', '82W6.068', '9N16.579',
  'Sebastian', '82W7.631', '9N15.274',
  'HS', '82W6.929', '9N16.773',
  'CK14', '82W7.557', '9N15.238',
  'CK4', '82W7.625', '9N15.517',
  'Moguls', '', ''
) %>%
  rowwise %>% 
  mutate(lat = convert_to_decimal(lat_garmin),
         lon = convert_to_decimal(lon_garmin)) %>%
  ungroup %>%
  select(-ends_with('garmin')) %>%
  # inner_join(site_list, 
  #            by = c('reef' = 'site')) %>% 
  filter(!is.na(lat)) %>%
  st_as_sf(coords = c('lon', 'lat')) %>%
  st_set_crs('WGS84')

#### Map ####
zoom_map <- bocas_shape %>%
  filter(!str_detect(folder, 'isla_colon')) %>%
  filter(!str_detect(basename, '[Rr]iver')) %>%
  filter(str_detect(basename, '[Ii]sla|Cayo|Region')) %>%
  # filter(str_detect(basename, 'Coast')) %>%
  ggplot(data = .) +
  geom_sf(fill= "antiquewhite") +
  # geom_sf(data = filter(tmp, str_detect(basename, 'Bath')) %>%
  #           mutate(depth = str_extract(basename, '[0-9]+') %>% as.integer,
  #                  depth_c = fct_reorder(as.character(depth), depth)),
  #         aes(colour = depth_c), alpha = 0.5) +
  geom_star(data = tibble(lon = -82.2408, lat = 9.3408),
            aes(x = lon, y = lat),
            fill = 'black', size = 4) +
  geom_sf(data = panama_locations) +
  scale_colour_brewer(palette = 'Greys', direction = -1) +
  coord_sf(xlim = c(-82.07, -82.25), ylim = c(9.2, 9.4), 
           expand = TRUE) +
  # coord_sf(xlim = c(-82.085, -82.22), ylim = c(9.25, 9.36), 
  #          expand = TRUE) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "aliceblue"))




carib_plot <- ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  geom_sf(data = tibble(lat = c(9.2, 9.4),
                        lon = c(-82.05, -82.25)) %>%
            st_as_sf(coords = c('lon', 'lat')) %>%
            st_set_crs('WGS84') %>%
            st_bbox() %>% 
            st_as_sfc() %>%
            st_as_sf(),
          fill = 'black') +
  # coord_sf(xlim = c(-88, -80), ylim = c(7.5, 10), 
  #          expand = TRUE) +
  coord_sf(xlim = c(-88, -60), ylim = c(7.5, 23),
           expand = TRUE) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "aliceblue"))



#### All in 1 ####
bocas_shape %>%
  filter(!str_detect(folder, 'isla_colon')) %>%
  filter(!str_detect(basename, '[Rr]iver')) %>%
  filter(str_detect(basename, '[Ii]sla|Cayo|Region')) %>%
  ggplot(data = .) +
  geom_sf(data = filter(world, name != 'Panama'), 
          fill= "antiquewhite") +
  geom_sf(fill= "antiquewhite") +
  geom_sf(data = panama_locations, size = 2) +
  coord_sf(xlim = c(-88, -60), ylim = c(7.5, 23),
           expand = TRUE) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  
  geom_magnify(from = c(xmin = -82.22, xmax = -82.085, ymin = 9.25, ymax = 9.36), 
               to = c(xmin = -82, xmax = -75, ymin = 12, ymax = 19),
               colour = 'black', shape = "rect",
               linewidth = 0.5, linetype = 'solid',
               proj.linetype = 'dotted', proj = 'facing',
               aspect = "fixed", expand = 0,
               corners = 0.25, shadow = FALSE) +
  
  annotation_scale(location = 'br', width_hint = 0.2, bar_cols = c('black', 'white'), 
                   line_col = 'black', text_col = 'black') + 
  annotation_north_arrow(location = 'tr', which_north = 'true', 
                         #pad_x = unit(0.75, 'in'), pad_y = unit(0.5, 'in'), 
                         style = north_arrow_fancy_orienteering(fill = c('grey25', 'gray75'), 
                                                                line_col = 'black', 
                                                                text_col = 'black')) 

#### Make as 2 and then insert for more complete control ####
#https://hughjonesd.github.io/ggmagnify/