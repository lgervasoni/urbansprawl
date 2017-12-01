###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################


####################################################################################
# Parameters
# Use a selection of leisure activities
LEISURE_ACTIVITIES_SELECTION = True

ACTIVITIES_CLASSIFICATION = ['shop', 'leisure/amenity', 'commercial/industrial']
####################################################################################

# Columns of interest corresponding to OSM keys 
columns_osm_tag = [ "amenity", "building", "landuse", "leisure", "shop" ]


'''
Classiy uses according to OpenStreetMap wiki
'''
####################################################################################
### Classification according to OSM wiki

### Key: amenity -> Related to activities (filtered)
amenities_sustenance = ['bar','pub','restaurant','biergarten','cafe','fast_food','food_court','ice_cream','pub','restaurant']
amenities_education = ['college','kindergarten','library','public_bookcase','school','music_school','driving_school','language_school','university']
amenities_transportation = ['fuel','bicycle_rental','bus_station','car_rental','taxi','car_wash','ferry_terminal']
amenities_financial = ['atm','bank','bureau_de_change']
amenities_healthcare = ['baby_hatch','clinic','dentist','doctors','hospital','nursing_home','pharmacy','social_facility','veterinary']
amenities_entertainment = ['arts_centre','brothel','casino','cinema','community_centre','fountain','gambling','nightclub','planetarium','social_centre','stripclub','studio','swingerclub','theatre']
amenities_others = ['animal_boarding','animal_shelter','courthouse','coworking_space','crematorium','dive_centre','dojo','embassy','fire_station','gym','internet_cafe','marketplace','police','post_office','townhall']

amenities_activities = amenities_sustenance + amenities_education + amenities_transportation + amenities_financial + amenities_healthcare + amenities_entertainment + amenities_others

### Key: shop -> all of them
shop_other = ["bookmaker","copyshop","dry_cleaning","e-cigarette","funeral_directors","laundry","money_lender","pawnbroker","pet","pyrotechnics","religion","tobacco","toys","travel_agency","vacant","weapons","user defined"]
shop_gifts = ["anime","books","gift","lottery","newsagent","stationery","ticket"]
shop_art = ["art","collector","craft","frame","games","model","music","musical_instrument","photo","camera","trophy","video","video_games"]
shop_sports = ["bicycle","car","car_repair","car_parts","fuel","fishing","free_flying","hunting","motorcycle","outdoor","scuba_diving","sports","swimming_pool","tyres"]
shop_electronics = ["computer","electronics","hifi","mobile_phone","radiotechnics","vacuum_cleaner"]
shop_furniture = ["antiques","bed","candles","carpet","curtain","furniture","interior_decoration","kitchen","lamps","tiles","window_blind"]
shop_household = ["agrarian","bathroom_furnishing","doityourself","electrical","energy","fireplace","florist","garden_centre","garden_furniture","gas","glaziery","hardware","houseware","locksmith","paint","security","trade"]
shop_health = ["beauty","chemist","cosmetics","drugstore","erotic","hairdresser","hairdresser_supply","hearing_aids","herbalist","massage","medical_supply","nutrition_supplements","optician","perfumery","tattoo"]
shop_charity = ["charity","second_hand","variety_store"]
shop_clothing = ["baby_goods","bag","boutique","clothes","fabric","fashion","jewelry","leather","shoes","tailor","watches"]
shop_mall = ["department_store","general","kiosk","mall","supermarket"]
shop_food = ["alcohol","bakery","beverages","brewing_supplies","butcher","cheese","chocolate","coffee","confectionery","convenience","deli","dairy","farm","greengrocer","ice_cream","organic","pasta","pastry","seafood","spices","tea","wine"]

shop_activities = shop_other + shop_gifts + shop_art + shop_sports + shop_electronics + shop_furniture + shop_household + shop_health + shop_charity + shop_clothing + shop_mall + shop_food + ['shop']

### Key: leisure 
if (not LEISURE_ACTIVITIES_SELECTION): # All leisure items
	leisure_activies = ['adult_gaming_centre','amusement_arcade','beach_resort','dance','hackerspace','ice_rink','pitch','sports_centre','stadium','summer_camp','swimming_area','water_park',"dog_park", "bird_hide", "bandstand", "firepit", "fishing", "garden", "golf_course", "marina", "miniature_golf", "nature_reserve", "park", "playground", "slipway", "track", "wildlife_hide", "swimming_pool"]
else: # Selection
	#Not tagged as activity: dog_park, bird_hide, bandstand, firepit, fishing, garden, golf_course, marina, miniature_golf, nature_reserve, park, playground, slipway, track, wildlife_hide, 
	leisure_activies = ['adult_gaming_centre','amusement_arcade','beach_resort','dance','hackerspace','ice_rink','pitch','sports_centre','stadium','summer_camp','swimming_area','water_park','swimming_pool']


### Key: building
building_other = ['barn','bridge','bunker','cabin','construction','cowshed','digester','garage','farm_auxiliary','greenhouse','hangar','hut','roof','shed','stable','sty','transformer_tower','service','ruins','yes','kiosk']
building_other_related_activities = ['kiosk', 'garages', 'hangar', 'stable', 'cowshed', 'digester'] + ['shop'] # From building_other related to activities
building_commercial = ['commercial','office','industrial','retail','warehouse'] + ['port']
building_civic_amenity = ['cathedral','chapel','church','mosque','temple','synagogue','shrine','civic','hospital','school','stadium','train_station','transportation','university','public']

building_activities = building_commercial + building_civic_amenity + building_other_related_activities
building_residential = ['hotel','farm','apartment','apartments','dormitory','house','residential','retirement_home','terrace','houseboat','bungalow','static_caravan','detached']

### Key: landuse
# Land usage: Activity
#landuse_activities = ['commercial','industrial','retail','port'] + ['quarry','salt_pond','construction','military','garages'] + ['shop']
landuse_activities = building_activities + shop_activities + amenities_activities + leisure_activies + ['quarry','salt_pond','construction','military']

# Land usage: Residential
landuse_residential = ['residential']

# Land usage not related to residential or activity uses
landuse_other = ['cemetery', 'landfill', 'railway']
landuse_water = ['water', 'reservoir', 'basin']
landuse_green = ['allotments','conservation', 'farmland', 'farmyard','forest','grass', 'greenfield', 'greenhouse_horticulture','meadow','orchard','pasture','peat_cutting','plant_nursery','recreation_ground','village_green','vineyard']
landuse_notResidentialActivity = landuse_other + landuse_water + landuse_green
####################################################################################