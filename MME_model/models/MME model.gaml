/**
 *  Fire Spread Model for Mediterranean and Middle Eastern Cities
 *  Authors: Yonatan Shaham (yjonas83@gmail.com), Itzhak Benenson (bennya@post.tau.ac.il)
 *  Licence: GPLv3.0 see: http://www.gnu.org/licenses/gpl-3.0.txt 
 */

model MME

global torus:false {
	//global
	int ignitions<-1;
	float step <- 1 #mn;
	bool create_trees<-true;
	float sim_start;
	float time_limit<-1 #hour;
	
	//meteorology
	int wind<-45 min:0 max:359; //wind direction 0=to the east
	list<float> wind_effect<-[1,0.75,0.5,0.3,0.2,0.3,0.5,0.75,1,0.75,0.5,0.3,0.2,0.3,0.5,0.75,1,0.75,0.5,0.3,0.2,0.3,0.5,0.75];
	list<int> wind_sectors;
	float bonus<-1.1;
	float wind_speed<-8.3 min:0 max:1200 #m/#s; 
	float pa<-1.205; //air density for brands at 20c
	float cd<-1.005; //air heat capacity for brands at 20c
	
	//burning
	int normal_temperature<-293; //20c = 20K+273K
	int max_temperature<-normal_temperature+1000;
	float brand_prob<-0.1 min:0 max:1; 
	int draft_sens<-65;  //sensetivity to wind for draft calculations, in degrees
	float burning_speed<-1;
	float impigement_range<-1;
	string veg_type<-"shrub";// among:["shrub","shrub-dry","trees_5m","trees_10m"];
	float veg_front_burn_time<-300#s;// among:[300#s,240#s,180#s,120#s];
	
	float L_general<-16; //mean fuel load default
	float trees_burn_range<-75 #m;
	float brand_veg_ignition_factor<-1;
	geometry trees;
	geometry impingemnt_area;
	geometry veg_radiation_area_60s;
	geometry veg_radiation_area_120s;
	geometry veg_radiation_area_180s;
	geometry veg_radiation_area_300s;
	geometry burning_area;
	float tree_top<-5; //height to trees
	geometry burning_windows;
	float trees_buffer<-4;  //allows to create buffer of vegetation around windows
	
	bool tuned<-true;
	float tune_factor<-1.7;
	
	list<room> rooms_w_windows;
	
	//GIS
	string CRS<-"EPSG:2039";
	string buildings_file<-"../includes/city/city_buildings.shp";
	file shape_buildings_file<-shape_file(buildings_file,CRS);
	string trees_file<-"../includes/city/city_trees_fine.shp";
	file shape_trees_file<-shape_file(trees_file,CRS);
	string bounds_file<-"../includes/city/city_bounds.shp";
	file shape_bounds_file<-shape_file(bounds_file,CRS);
	geometry shape <- envelope(shape_bounds_file);	
	
	string rooms_file<-"../includes/city/corridor.shp";
	file shape_rooms_file<-shape_file(rooms_file,CRS);
	string doors_file<-"../includes/city/city_doors.shp";
	file shape_doors_file<-shape_file(doors_file,CRS);
	string windows_file<-"../includes/city/city_windows.shp";
	file shape_windows_file<-shape_file(windows_file,CRS);
	
	float balcony_ratio; // portion of windows with balconies
	
	list<building_line> locations;
	float floor_height<-2.8 #m;
		
	//STATISTICS AND OUTPUTS
	bool report_init<-true;
	
	string output_file_name<-"../doc/output.csv";
		
	string origin_of_fire<-nil;
	int origin_of_fire_building<-nil;
	int origin_of_fire_floor<-nil;
	point point_of_origin;
	
	//used to recored evey 5 minutes
	list<float> burnt_area_series<-nil;
	list<int> burning_buildings_series;
	list<int> burning_floor_series;
	list<int> burning_apt_series;
	list<int> burning_room_series;
	
	//fire spread mechanisms documentation
	bool brand_were_released<-false;
	float time_to_brands;
	bool veg_was_ignited<-false;
	float time_to_veg;
	bool ignited_by_rad<-false;
	float time_to_rad;
	point veg_ignition;
	int top_burning_floor<-0;
	
	//spread documentation
	float time_to_first_FO;
	float time_to_next_apt;
	float d_to_next_apt;
	room next_ignited_apt;
	float time_to_next_floor;
	float d_to_next_floor;
	building_line next_ignited_floor;
	float time_to_next_building;
	float d_to_next_building;
	int angel_to_next_building;
	string next_ignited_building;
	bool polygons_recording<-false;
	bool save_polygones<-false;
	
	//serial experiment
	int to_explore<-0;
	list<string> apts; 
	int serial<-0;
	list<int> UNIQS;
	string room_name;
	int floor_number;
	string output_origin_file_name;
	bool go_serial;
	bool serial_building;
	float step_time;
	float zero_time;
	
	init{
		float init_time<-machine_time;
		float start_time<-init_time;
		zero_time<-init_time;
		write "Starting init phase.";	
		create room from: shape_rooms_file with: [floors::int(read("floors")),floor::1,UNIQ_ID::int(read("UNIQ_ID")),apt::int(read("apt")),Ltag::float(read("L"))]; 
		create door from: shape_doors_file with:[floor::1];
		create window from: shape_windows_file with:[floor::1,width::float(read("width")),height::float(read("height")),heading::float(read("heading"))];
		create building_line from:shape_buildings_file  with:[floor::int(1),UNIQ_ID::int(read("UNIQ_ID")),flame_roof::int(read("roof"))];
		
		
		
		write "number of buildings: "+length(building_line);
		ask room{
			Arf<-shape.area ;
			L<-Ltag*Arf #kg; 
		} 
		
		if report_init {write "GIS read took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		ask building_line{
			my_rooms<-room overlapping (self+0.1);
			floors<-max(my_rooms accumulate each.floors);
		}
		
		ask window {Ad<-height*width; location<-{location.x,location.y,0};}
		ask door{
			list<room> temp_rooms<-room overlapping (self+0.1);
			int max_floors<- max(temp_rooms accumulate each.floors);
			loop i from:2 to: max_floors{
				create door number:1 {
					floor<-i;
					shape<-myself.shape.contour;
					shape<-point(remove_duplicates(shape.points));
					location<-myself.location;
				}
			}
		}
		if report_init {write "doors took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;	
		
		ask window{
			list<room> temp_rooms<-room overlapping (self+0.1);
			int max_floors<- max(temp_rooms accumulate each.floors);
			
			loop i from:2 to: max_floors{
				create window number:1 {
					floor<-i;
					shape<-myself.shape.contour;
					shape<-point(remove_duplicates(shape.points));
					location<-myself.location;
					height<-myself.height;
					width<-myself.width;
					heading<-myself.heading;
					Ad<-height*width;
				}
			}
		}
		if report_init {write "windows took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		ask room{
			
			building_line my_building <-first(building_line where (each.floor=1 and each.UNIQ_ID = UNIQ_ID));
			if floors>1{
				loop i from:2 to:floors {
					create room number:1 {
						child<-true;
						floor<-i;
						floors<-myself.floors;
						shape<-myself.shape.contour;
						shape<-polygon(remove_duplicates(shape.points));	
						apt<-myself.apt;
						
						
						
						location<-myself.location;
						UNIQ_ID<-myself.UNIQ_ID;
						wr<-shape.height;
						dr<-shape.width ;
						Arf<-shape.area ; 	
						L<-Ltag*Arf ; 
						
						float calc_area<-wr*dr;
						float correction_ration<-Arf/calc_area;
						dr<-dr*correction_ration^0.5;
						wr<-wr*correction_ration^0.5;
											
						if i=myself.floors and my_building.flame_roof=1 {
							flame_roof<-true;
						}
						
					}
				}
			} else {
				if my_building.flame_roof=1 {
					flame_roof<-true;
				}
			}
			
			
			if floor=1 and my_building.floors=0 and my_building.flame_roof=1 {
				flame_roof<-true;
			}
		}
		
	
		if report_init {write "rooms took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		ask building_line{
			if floors>1{
				loop i from:2 to:floors{
					create building_line number:1{
						floor<-i;
						floors<-myself.floors;
						flame_roof<-myself.flame_roof;
						shape<-myself.shape.contour;
						shape<-polygon(remove_duplicates(shape.points));	
						location<-myself.location;
						UNIQ_ID<-myself.UNIQ_ID;
					}			
				}
			}
		}
		
		if report_init {write "floors took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		//create blaconies
			list<window> temp<-window;
			temp<-shuffle(temp);
			temp<-copy_between(temp,0,length(temp)*balcony_ratio);
			ask temp {
				balcony<-true;
			}
		//end balconies
		if report_init {write "balconies took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		loop i over: window{
			list<window> temp<-nil;
			ask i{
				 temp<-((window at_distance 5) where ((each.floor=floor) and each.location overlaps (self+0.01))) - self;
				}
			if length(temp)>0 {ask temp {do die;}}
		}
		if report_init {write "SPAT windows1 took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		ask room {
			float search_d<-max([shape.width*1.1/2,shape.height*1.1/2]);
			my_windows<-((window at_distance search_d) where ((each.floor=floor) and each.location overlaps (self+0.01)));
			
			my_doors<-((door at_distance search_d) where (each.floor=floor)) overlapping(self+0.1);
			ask my_doors{
				my_rooms<-remove_duplicates(my_rooms+myself);
			}
			ask my_windows{
				my_room<-myself;
			}
			top<-floor_height;
		}
		if report_init {write "SPAT rooms1 took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		ask building_line {my_rooms<-((room at_distance 10) where (each.floor=floor)) overlapping (self+0.1);}
		if report_init {write "SPAT buildings took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		ask window {
			heading_check<-heading;
			internal_heading<-(one_of(my_room) towards self); //starting to caculate real heaing (GIS is only 0<90)
			if heading<0 {heading<-heading+360;}
			if heading=360 {heading<-0;} //end of caculating real heading
			if wind_speed>=(5) {
				if ((internal_heading+180-wind<draft_sens) and (internal_heading+180-wind>0)) or ((wind-internal_heading+180<draft_sens) and (wind-internal_heading+180>0)) {upwind<-true;} //determining realtion to wind 
				if (internal_heading>270) and ((internal_heading-180-wind<draft_sens) and (internal_heading-180-wind>0)) or ((wind-internal_heading-180<draft_sens) and (wind-internal_heading-180>0))  {upwind<-true;}
				if ((internal_heading-wind<draft_sens) and (internal_heading-wind>0)) or ((wind-internal_heading<draft_sens) and (wind-internal_heading>0))  and !upwind {downwind<-true;}
			}
		}
		
		if report_init {write "SPAT windows2 took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		ask room { //setting room draf conditions and fire parmeters
			if int(my_windows count each.upwind)>=1 and  int(my_windows count each.downwind)>=1 {draft<-true;}
			Ad<- sum(my_windows accumulate (each.Ad));
			my_hd<- mean(my_windows accumulate (each.height));
			if my_hd=0 {my_hd<-0.762;} //half a door
			if Ad=0 {Ad<-0.762*0.32/2;}
			ArT<-(wr*2+dr*2)*floor_height+wr*dr*2-Ad;
			psy<-L/(Ad*ArT)^0.5;
			eta<-ArT/(Ad*my_hd^0.5);
			
			if draft{
				mr<-L/1200;
				
				Lt30<-L/mr/60*0.3;
				Lt80<-L/mr/60*0.8;
				Lt_total<-L/mr/60;
				
				
				Tr<-normal_temperature+1200*(1-#e^(-0.04*psy));
			}
			
			if !draft{
				mr<-0.18*(1-#e^(-0.036*eta))*Ad*(my_hd^0.5)*((dr/wr)^(-0.5));
				if tuned and length(my_windows)>=2{
					mr<-mr*tune_factor;
				}
				
					Lt30<-L/mr/60*0.3;
					Lt80<-L/mr/60*0.8;
				
				
				Lt_total<-L/mr/60;
				Tr<-normal_temperature+6000*(1-#e^(-0.1*eta))*(1-#e^(-0.05*psy))*(eta^(-0.5));
			}
		}
		
		if report_init {write "draft took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		
		
		//TREES
		if create_trees {
			create tree from:shape_trees_file with:[top::tree_top] {
				trees<-trees+shape;
				impingemnt_area<-impingemnt_area+(shape+impigement_range);
				
			}
		}	
		if trees_buffer>0{
			trees<-trees - (geometry(window) + trees_buffer);
			impingemnt_area<-impingemnt_area - (geometry(window)+ trees_buffer-impigement_range);
			
		}
		
		
		
		if report_init {write "trees took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		//MOVING FLOORS UP
		ask room{location<-location+{0,0,(floor-1)*floor_height};} //putting rooms in floors
		ask door{location<-location+{0,0,(floor-1)*floor_height};}
		ask window{
			location<-location+{0,0,(floor-0.5)*floor_height-height/2};
		}
		ask building_line {location<-location+{0,0,(floor-1)*floor_height};}
		
		ask room {apt_id<-string(UNIQ_ID)+"_"+string(floor)+"_"+string(apt);}
		ask building_line {floor_id<-string(UNIQ_ID)+"_"+string(floor);}
		
		//killing erroed doors and windows
		ask window where (length(each.my_room)=0) {do die;}
		ask door where (length(each.my_rooms)=0) {do die;}
		
		if report_init {write "remove agents took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		rooms_w_windows<-room where (each.L>0 and length(each.my_windows)>0 and each.floor<=2) sort_by each.name;
		locations<-building_line where (each.floor=1);
		
		apts<-remove_duplicates((room where (each.floor=1 and each.apt>0)) accumulate each.apt_id);
		write "number of apartments: "+length(apts);
		if report_init {write "floors UP took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
		
		//radiation line of sight clculations
		ask window{
			in_sight<-(window at_distance 10) where (empty((room) overlapping line([location,each.location])));
			in_sight<-in_sight where (first(my_room).UNIQ_ID!=first(each.my_room).UNIQ_ID);
					
		}
		if report_init {write "windwos in sight took :" + ((machine_time-init_time)/1000);}
		init_time<-machine_time;
				
		ask room {UNIQ_ID<-UNIQ_ID;}
		UNIQS<-remove_duplicates(room accumulate each.UNIQ_ID) sort_by each;
		ask window { UNIQ_ID<-first(my_room).UNIQ_ID; }
		if report_init { write "init phase took:" + ((machine_time-start_time)/1000); write "=====";}
		init_time<-machine_time;
		
	}
	
	action project_burnt_area (float craetion_timing){
		geometry temp<-nil;
		loop i over:building_line where (each.color=#red or each.color=#black){
			building_line temp_room<-i;
			loop j from: 0 to: length(temp_room.shape.points) - 1{ temp_room.shape <-  set_z (temp_room.shape, j, 0);}
			temp<-temp union temp_room.shape;
		}
		temp<-temp union burning_area;
		geometry older_projections<-union(projected_burn_area accumulate each.shape);
		create projected_burn_area number:1 {shape<-temp; timing<-craetion_timing;}
	}
	
	reflex go {
		sim_start<-machine_time;
		step_time<-machine_time;
		//Wind roses from (Albini, 1976)
		if wind_speed=4{
			wind_effect <- [1,0.57,0.28,0.18,0.16,0.18,0.28,0.57,1,0.57,0.28,0.18,0.16,0.18,0.28,0.57,1,0.57,0.28,0.18,0.16,0.18,0.28,0.57];
			
		}
		
		if wind_speed=6{
			wind_effect <- [1,0.49,0.25,0.14,0.13,0.14,0.21,0.49,1,0.49,0.25,0.14,0.13,0.14,0.21,0.49,1,0.49,0.25,0.14,0.13,0.14,0.21,0.49];
			
		}
		
		if wind_speed=8{
			wind_effect <- [1,0.25,0.09,0.04,0.04,0.04,0.09,0.25,1,0.25,0.09,0.04,0.04,0.04,0.09,0.25,1,0.25,0.09,0.04,0.04,0.04,0.09,0.25];
			
		}
		
		if veg_type="shrub"{
			switch wind_speed{
				match 4.0 {
					burning_speed <- 0.011;
				}
				match 6.0 {
					burning_speed <- 0.019;
				}
				match 8.0 {
					burning_speed <- 0.027;
				}
			}
		}
		
		if veg_type="shrub-dry"{
			switch wind_speed{
				match 4.0 {
					burning_speed <- 0.016;
				}
				match 6.0 {
					burning_speed <- 0.027;
				}
				match 8.0 {
					burning_speed <- 0.052;
				}
			}
		}
		if time=0{
			if go_serial{ //setting up ignition in ground floor apartments one by one
				list<room> pot_rooms<-room where (each.apt_id=apts[to_explore] and each.floor=1 and length(each.my_windows)>0 and each.L>0);
				ask one_of(pot_rooms){
					do ignite;
					origin_of_fire<-apt_id;
					origin_of_fire_building<-UNIQ_ID;
					origin_of_fire_floor<-floor;
					point_of_origin<-{location.x,location.y,0};
					time_to_first_FO<-Lt30;
					}	
			} else {
				if serial_building{
					list<room> pot_rooms<-room where (each.UNIQ_ID=UNIQS[to_explore] and  length(each.my_windows)>0 and each.L>0);
					
					ask one_of(pot_rooms){
						do ignite;
						origin_of_fire<-apt_id;
						origin_of_fire_building<-UNIQ_ID;
						origin_of_fire_floor<-floor;
						point_of_origin<-{location.x,location.y,0};
						time_to_first_FO<-Lt30;
						}		
				}
					
			}
		}
		
		
		//trees burning
		ask seg {do die;}
		
		wind_sectors <- nil;
		loop i from:0 to:8*3{
			int temp <-wind -360 + 45*(i) - 22;
			
			wind_sectors <- wind_sectors + temp;
		}
		
		if burning_area!=nil{
			list<point> seg_points <- points_on(burning_area.contour,1);
			if length(seg_points)>10{
				loop i from:0 to: length(seg_points)-2{
				create seg number:1 {
					shape <- line([seg_points[i],seg_points[i+1]]);
				}
				}
				create seg number:1 {
					shape <- line([last(seg_points),first(seg_points)]);
				}
			}
		}
		
		//
		ask seg{do set_perpendicular;}
		list<point> targets <- seg collect each.target;
		if length(targets)>8 {
			
			burning_area <- polygon(targets);
		} else {
			burning_area<-burning_area+burning_speed*step;
		}
		burning_area<-burning_area intersection impingemnt_area;
		burning_windows<-nil;
		if report_init {write "trees took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;
		
		//trees radiation
		
		if (veg_type="shrub" or veg_type="shrub-dry") and veg_front_burn_time=300 {
			veg_radiation_area_300s<-burning_area+3-geometry(building_line where (each.floor=1));
			veg_radiation_area_300s<-veg_radiation_area_300s+0.1;
		}
		if veg_type="trees_5m"{
			switch veg_front_burn_time{
				match 120.0 {
					veg_radiation_area_180s<-burning_area+7-geometry(building_line where (each.floor=1));
					veg_radiation_area_120s<-burning_area+6-geometry(building_line where (each.floor=1));
					veg_radiation_area_60s<-burning_area+5-geometry(building_line where (each.floor=1));
					
					veg_radiation_area_180s<-veg_radiation_area_180s+0.1;
					veg_radiation_area_120s<-veg_radiation_area_120s+0.1;
					veg_radiation_area_60s<-veg_radiation_area_60s+0.1;
				}
			}
		}	
		
		
		//trees ignition of windows
		list<window> toched_by_fire <- window where (({each.location.x,each.location.y} overlaps burning_area) and (each.location.z<=tree_top*(1+1/3)));
		ask toched_by_fire {ask my_room {do ignite;}}
		list<window> not_effected_by_veg_rad <-window where (each.ignition_time_by_veg_rad=0 and each.location.z<=10);
	
		ask not_effected_by_veg_rad where (({each.location.x,each.location.y} overlaps veg_radiation_area_300s))  {do recive_rad_from_veg(300);}
		ask not_effected_by_veg_rad where (({each.location.x,each.location.y} overlaps veg_radiation_area_180s)) {do recive_rad_from_veg(180);}
	
		ask not_effected_by_veg_rad where (({each.location.x,each.location.y} overlaps veg_radiation_area_120s)) {do recive_rad_from_veg(120);}
		ask not_effected_by_veg_rad where (({each.location.x,each.location.y} overlaps veg_radiation_area_60s)) {do recive_rad_from_veg(60);}
		
		ask window where (each.ignition_time_by_veg_rad>0){
			
			if time>=ignition_time_by_veg_rad{
				ask my_room[0] {do ignite;}
				ignition_time_by_veg_rad<--1;
				
			}
		}
		if report_init {write "windows took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;
		
		
		
		if report_init {write "NON-BURNING rooms took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;
		ask (room where each.burn) {do update;} 
		if report_init {write "BURNING rooms took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;

		//fire spread by brands
		ask brand{
			do ignite;
			do die;
		}
		if report_init {write "brands took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;
		
		ask building_line{
			if length(my_rooms where ((each.burn) and !(each.burnt_down)))>0 {color<-#red;}
			if length(my_rooms where (!(each.burn) and (each.burnt_down)))>0 and length(my_rooms where ((each.burn) and !(each.burnt_down)))=0 {color<-#black;}
		}
		if report_init {write "building lines took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;
		ask roof_flame {do update;}
		if report_init {write "roof flames took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;

		//STATISTICS AND OUTPUT		
		int max_buildings<-max(burning_buildings_series);
		int max_floors<-max(burning_floor_series);
		int max_apts<-max(burning_apt_series);
		burnt_area_series<-burnt_area_series+sum(((room where ((each.burn) and !(each.burnt_down))) + (room where (!(each.burn) and (each.burnt_down)))) accumulate each.shape.area);
		burning_buildings_series<-burning_buildings_series + length(remove_duplicates((building_line where ((each.color=#red))) accumulate each.UNIQ_ID));
		burning_floor_series<-burning_floor_series + length(remove_duplicates((building_line where ((each.color=#red))) accumulate each.floor_id));
		burning_apt_series<-burning_apt_series + length(remove_duplicates((room where ((each.burn) and !(each.burnt_down))) accumulate each.apt_id));
		burning_room_series<-burning_room_series + length(room where ((each.burn) and !(each.burnt_down)));
		
		if !(burning_area=nil) {
			if burning_area.area>0 {
				if !veg_was_ignited {time_to_veg<-time;}
				veg_was_ignited<-true;
			}
		}
		if max_buildings=ignitions and max(burning_buildings_series)>ignitions {
			time_to_next_building<-time;
			
			room ignited <- room where ((each.burn) and each.UNIQ_ID!=origin_of_fire_building) closest_to point_of_origin;
			write "ignited " + ignited;
			next_ignited_building <- ignited.UNIQ_ID;
			point new_building <- ignited;
			d_to_next_building <- point_of_origin distance_to new_building;
			angel_to_next_building <-  new_building towards point_of_origin  ;
			write "First building ignited! origin: " + origin_of_fire_building + ", target: " + next_ignited_building + ", d: " + d_to_next_building + ", angel_to_next_building: "+angel_to_next_building; 
		}
		if max_floors=ignitions and max(burning_floor_series)>ignitions {
			time_to_next_floor<-time;
		}
		if max_apts=ignitions and max(burning_apt_series)>ignitions {
			time_to_next_apt<-time;
		}
		top_burning_floor<-max((building_line where ((each.color=#red))) accumulate each.floor);
		if report_init {write "STATS took :" + ((machine_time-sim_start)/1000);}
		sim_start<-machine_time;
		if report_init {write "ALL took :" + ((machine_time-step_time)/1000);}
		
		if report_init {write "----------";}
		
		if time>0 { //STOPPING simulation	
			if time>=time_limit{
				
				do pause;
			} 
		}
	}
}

species seg{
	int secotr<--1;
	point target;
	float factor;
	int angle;
	int perpendicular;
	
	action set_perpendicular{
		angle <- first(shape.points) towards last(shape.points);
		perpendicular <- angle+90;
		if perpendicular > 360 {perpendicular<-perpendicular-360;}
		//set sector factor
		float spread_length<-burning_speed*step;
		int sector;
		float sector_factor<-0;
		loop i from:0 to:15{
			if between(perpendicular,wind_sectors[i],wind_sectors[i+1])  {
				
				sector_factor <- wind_effect[i];
				sector<-i;
			}
		}
		if perpendicular<first(wind_sectors) and perpendicular>last(wind_sectors) {
			write ""+self+" overloop";
			sector_factor <- wind_effect[0];
				sector<-0;
		}
		spread_length <- spread_length * sector_factor * bonus;
		target <- self translated_by {spread_length*cos(perpendicular),spread_length*sin(perpendicular)};
	}
	
	
	
	aspect base{
		draw shape color:#blue;
		draw polygon([first(shape.points),target,last(shape.points)]) color:#red;
	}
}

species roof_flame{
	int UNIQ_ID;
	list<room> my_rooms;
	int floor;
	
	float Af;
	float Pr<-0.27;
	
	float Df;
	float Hf;
	float f_teta;
	
	float u_ast;
	float mrf;
	float yr;
	float wd;
	float e_f;
	
	list<building_line> radiated_roofs;
	list<building_line> in_range;
	list<window> radiated_windows;
	
	action update{
		my_rooms<-room where (each.UNIQ_ID=UNIQ_ID and (each.floor=each.floors or each.floor=1) and each.flame_roof and each.burn);
		if !empty(my_rooms){
			Af<-sum(my_rooms accumulate each.Arf);
			Df<-(4*Af/#pi)^0.5/2;
			yr<-mean(my_rooms accumulate (floor_height/2-each.my_hd));
			wd<-mean(my_rooms accumulate (each.Ad/each.my_hd));
			mrf<-0.383*yr^1.5*wd;
			u_ast<-wind_speed/(9.8*mrf*Df/Pr)^(1/3);
			Hf<-Df*55*((mrf/Pr/(9.8*Df)^0.5)^0.67)*(u_ast)^(-0.21);
			f_teta<- u_ast<1 ? acos(1) : acos(1/(u_ast^0.5));
			e_f<-140*exp(-0.12*Df)+20*(1-exp(-0.12*Df));
			
			if wind_speed<5{
				radiated_roofs<-(building_line at_distance 100) where (each.floor>=floor and each.flame_roof and each.floors=each.floor);
			} else {
				radiated_roofs<-(building_line at_distance 100) where (each.floor>=floor and each.flame_roof and each.floors=each.floor and ((each towards self)>=wind+90 or (each towards self)<=wind-90));
			}
		
			ask radiated_roofs {
				float L<-{location.x,location.y,0} distance_to {myself.location.x,myself.location.y,0};
				
				if abs(myself.location.z-location.z)<1{
					float S<-2*L/myself.Df;
					float h<-2*myself.Hf/myself.Df;
					float A<-(h^2+S^2+1)/(2*S);
					float B<-(1+S^2)/(2*S);
					float F<-(B-1/S)/(#pi*(B^2-1)^0.5)*atan(( (B+1)*(S-1)/(B-1)*(S+1)  )^0.5)  - (A-1/S)/(#pi*(A^2-1)^0.5)*atan(( (A+1)*(S-1)/(A-1)*(S+1)  )^0.5);
					q<-myself.e_f*F;
					do recive_radiation(q);
				} else {
					if myself.location.z-location.z<0 {
						float h1<-location.z-myself.location.z;
						float h2<-myself.Hf-h1;
						float S<-2*L/myself.Df;
						float h<-2*(h2)/myself.Df;
						float A<-(h^2+S^2+1)/(2*S);
						float B<-(1+S^2)/(2*S);
						float F<-(B-1/S)/(#pi*(B^2-1)^0.5)*atan(( (B+1)*(S-1)/(B-1)*(S+1)  )^0.5)  - (A-1/S)/(#pi*(A^2-1)^0.5)*atan(( (A+1)*(S-1)/(A-1)*(S+1)  )^0.5);
						q<-myself.e_f*F;
						do recive_radiation(q);
					}
				}
			}
		
			float z_point<-location.z;
			
			if wind_speed<5{
				in_range<-(building_line at_distance 100) where (each.floor>=floor);
			} else {
				in_range<-(building_line at_distance 100) where (each.floor>=floor and ((each towards self)>=wind+90 or (each towards self)<=wind-90));
			}
		
			ask in_range{
				if myself.location.z-location.z<(-1) {
					myself.radiated_windows<-(window at_distance 20) where (each.UNIQ_ID=UNIQ_ID and each.location.z>=z_point);
				}
			}
		
			radiated_windows<-radiated_windows where (empty((room + building_line + window) overlapping line([location+{0,0,1},each.location])));
			ask radiated_windows{
				float L<-{location.x,location.y,0} distance_to {myself.location.x,myself.location.y,0};
				
				if abs(myself.location.z-location.z)<1{
					float S<-2*L/myself.Df;
					float h<-2*myself.Hf/myself.Df;
					float A<-(h^2+S^2+1)/(2*S);
					float B<-(1+S^2)/(2*S);
					float F<-(1/#pi/S) * atan(h / (S^2-1)^0.5 ) - (h/S/#pi)*atan( ( (S-1)/(S+1) )^0.5 ) + (A*h/#pi/S/( A^2-1 )^0.5)*atan( ( (A+1)*(S-1)/(A-1)/(S+1) )^0.5 );
					float q<-myself.e_f*F;
					do recive_radiation(q);
				} else {
					if myself.location.z-location.z<0 {
						float h1<-location.z-myself.location.z;
						float h2<-myself.Hf-h1;
						float S<-2*L/myself.Df;
						float h<-2*(h2)/myself.Df;
						float A<-(h^2+S^2+1)/(2*S);
						float B<-(1+S^2)/(2*S);
						float F<-(1/#pi/S) * atan(h / (S^2-1)^0.5 ) - (h/S/#pi)*atan( ( (S-1)/(S+1) )^0.5 ) + (A*h/#pi/S/( A^2-1 )^0.5)*atan( ( (A+1)*(S-1)/(A-1)/(S+1) )^0.5 );
						float q<-myself.e_f*F;
						S<-2*L/myself.Df;
						h<-2*(h1)/myself.Df;
						A<-(h^2+S^2+1)/(2*S);
						B<-(1+S^2)/(2*S);
						F<-(1/#pi/S) * atan(h / (S^2-1)^0.5 ) - (h/S/#pi)*atan( ( (S-1)/(S+1) )^0.5 ) + (A*h/#pi/S/( A^2-1 )^0.5)*atan( ( (A+1)*(S-1)/(A-1)/(S+1) )^0.5 );
						do recive_radiation(q);
					}
				}
			}
		} else {
			do die;
		}
	}
	
	aspect base_3D{
		draw cylinder(Df,Hf) color:#yellow;
		draw line([self,{self.location.x+cos(wind)*cos(90-f_teta)*Hf,self.location.y+sin(wind)*cos(90-f_teta)*Hf,self.location.z+sin(90-f_teta)*Hf}]) color:#black;
	}
}

species projected_burn_area{
	float timing;
	
	aspect base{
		draw shape color:#purple empty:true border:#purple;
	}
}

species building_line{
	int UNIQ_ID;
	float BLDG_HT;
	int floors;
	int floor;
	float burning_area;
	float burnt_area;
	list<room> my_rooms;
	rgb color<-#white;
	string floor_id;
	int flame_roof;
	
	//radiation
	float q;
	float radiation_time;
	
	action recive_radiation(float influx){
		switch influx{
			match_between [0,12500] {radiation_time<-0;}
			match_between [12500,15000] {radiation_time<-radiation_time + step;}
			match_between [15000,17500] {radiation_time<-radiation_time + step*30/25;}
			match_between [17500,20000] {radiation_time<-radiation_time + step*30/10;}
			match_between [20000,30000] {radiation_time<-radiation_time + step*30/7;}
			match_between [30000,3000000] {radiation_time<-radiation_time + step*30/1;}
		}
		if radiation_time>30*60 {
			ask one_of(my_rooms where (length(each.my_windows)>1)) {do ignite;}
			ignited_by_rad<-true;
			time_to_rad<-time;
		}
	}
	
	aspect base_3D{
		draw shape color: color depth:floor_height;
		if flame_roof=1 and floor=1 {draw shape+2 color: #gold;}
	}
	aspect base{
		if flame_roof=1  and floor=1 {draw shape+2 color: #gold;}
		draw shape color: color ;
		if length((building_line at_distance 5) where (each.color=#red and each.UNIQ_ID=UNIQ_ID))>0 {
			draw shape ;
		}
		if length((building_line at_distance 5) where (each.color=#black and each.UNIQ_ID=UNIQ_ID))>0 {
			draw shape;
		}
	}
}

species element{
	rgb color <- #white;
	float buttom<-0.0;
	float top<-10.0;
	int floor;
	float ingition_time<-nil;
	bool burn<-false;
	bool burnt_down<-false;
	int temperature<-normal_temperature;
	int burning_for<-0;
	
	

	action change_to_red {color<-#red;}
	
	action set_color{
		color <-rgb(255,255*(max_temperature-temperature)/(max_temperature-normal_temperature),255*(max_temperature-temperature)/(max_temperature-normal_temperature));
		if burnt_down{color<-#black;}
	}
	
	action ignite  {
		if !burn and !burnt_down{
			burn<-true;
			do set_color;
			ingition_time<-time;
			} 
	}
	
	action update {
		
	}
	
	aspect base {
		draw shape color: color;			
	}
}

species room parent:element{
	list<list> my_rooms;
	list<window> my_windows<-nil;
	list<door> my_doors<-nil;
	bool ignited_by_door<-false;
	bool draft<-false;
	int number;
	int UNIQ_ID;
	float BLDG_HT<-nil;
	int apt;
	int floors;
	string apt_id;
	bool child;
	//Parameters for Lee and Davidson model of fire development (Lee & Davidson, 2010)

	float width;
	float height;
	float wr<-shape.height #m; 	//will be chacked for replacement during init
	float dr<-shape.width #m; 		//will be chacked for replacement during init
	float Ltag<-L_general #kg/#m2; 		//fuel load density
	float Arf<-shape.area #m2; 		// floor area
	float L<-Ltag*Arf #kg; 			//total fuel load
	float Ad<-nil;				// window area, calculaten in global init
	float ArT<-nil;				//surface area minus window
	float psy<-nil;
	float eta<-nil;
	float my_hd<-nil;
	bool flame_roof<-false; 	//is the roof flamable ?
	
	//variables for Lee and Davidson model of fire development (Lee & Davidson, 2010)
	float mr<-nil;				//rate of burning
	float Lt30<-nil;			//flash over start time (minutes)
	float Lt80<-nil;			//flash over end time (minutes)
	float Lt_total<-nil;		//total burn time (minutes)
	float Tr<-nil;				//flash over tempetature (K)
	float flash_start<-nil;		//flash over absolut start time (seconds)
	float flash_end<-nil;		//flash over absolut end time (seconds)
	float burning_end<-nil;		//total absolut burn time (seconds)
	float brands_number<-nil;
	
	
	
	action ignite  {
		if !burn and !burnt_down{
			burn<-true;
			do set_color;
			ingition_time<-time;
			flash_start<-ingition_time+Lt30*60;
			flash_end<-ingition_time+Lt80*60;
			burning_end<-ingition_time+Lt_total*60;
			} 
	}
	
	action update {
		
		if burn and time>=flash_start and time<flash_end { //develop fire section
			temperature<-Tr; 
			do change_to_red;
			if burning_for=0 {ask my_doors {do impige;}}
			ask my_windows {do impige; do radiate;}
			if flame_roof {
				brands_number<-Arf*306.77*exp(0.1879*wind_speed)*step/(flash_end-flash_start);
				
				loop times:round(brands_number*(0.27*0.005+0.02*0.02)){ //probability that brand is effective for roofs
					do release_barnds;
				}
				
				if empty(roof_flame where (each.UNIQ_ID=UNIQ_ID)){
					create roof_flame number:1 {
						location<-first(building_line where (each.UNIQ_ID=myself.UNIQ_ID and each.floor=myself.floor)).location;
						location<-location+{0,0,floor_height};
						my_rooms<-remove_duplicates(my_rooms+myself);
						UNIQ_ID<-myself.UNIQ_ID;
						floor<-myself.floor;
					}	
				} 
			}
		
			burning_for<-burning_for+1;
			} 
		
		if burn and ((time<flash_start) or (time>=flash_end and time<burning_end)){
			temperature<-Tr/2;
		}
		
		//write ""+self+
		do set_color;
		
		if burn and time>=burning_end {
			burn<-false;
			burnt_down<-true;
			color<-#black;
			temperature<-normal_temperature;
		}
	}
	
	action release_barnds{
			if !brand_were_released{
				time_to_brands<-time;
			}
			brand_were_released<-true;
			float uLx<-nil; //Himoto and Tanaka adopted by Lee and Davidson
			float SigLx<-nil;
			float B<-nil;
			float pb<-125;	
			float db<-1.4;	
			
			B<-wind_speed*(9.8*Arf^0.5)*((pb*db)/(pa*Arf^0.5))^(-3/4)*((1500*Arf)/(pa*cd*normal_temperature*(9.8^0.5)*Arf^(0.5*5/2)))^0.5;
			SigLx<-(ln(1+(0.88*B^(1/3)/0.47*B^(2/3)))^2)^0.5;
			uLx<-ln(Arf^0.5*(0.47*B^(2/3))/(1+(0.88*B^(1/3)+0.47*B^(2/3))^2)^0.5);
			
			float px<-(abs(gauss(uLx,SigLx)));
			float py<-gauss(0,0.92*Arf^0.5);
			float right<-nil; //wiil the brand fall right or left from main axis
			if flip(0.5){right<-1;} else {right<-(-1);}
			float Bteta<-nil;
			if px=0 {Bteta<-90*right+wind;} else {Bteta<-atan(py/px)*right+wind;}
			float Bdistance<-(px^2+py^2)^0.5;
			create brand {
				location<-myself.location+{cos(Bteta)*Bdistance,sin(Bteta)*Bdistance};
				location<-{location.x,location.y,0};
			}
	}
	
	user_command ignite_self {  //allows user to manually ignite the room
    	do ignite;
    	point_of_origin <- self;
	}
	
	aspect over_view{
		draw shape color: color depth:floor_height;
		if flame_roof {draw shape translated_by {0,0,floor_height+0.01} color:#gold;}
	}
	
	aspect base {
		draw shape color: color border:#black;
	}
}


species window parent:element{
	float heading<-nil; 
	list<room> my_room<-nil;
	int UNIQ_ID;
	float internal_heading<-nil;
	float heading_check<-nil;
	bool upwind<-false;
	bool downwind<-false;
	float width;
	float height;
	float Ad<-nil;
	float xw<-nil; 					//L&D model flame center line projection distance
	float flame_impigment_range;
	float flame_temperature<-813; 	//L&D model
	float flame_tip<-nil;
	bool balcony<-false;
	geometry flame<-nil;
	geometry casted_flame<-nil;
	
	float ignition_time_by_veg_rad;
	
	//radiation
	list<window> in_sight;
	bool show_in_sight;
	float temp_rad;
	float ratiation_level;
	float radiation_time;
	
	//imping and rad pre-calculations
	bool calculated;
	map<window,float> factors;
	map<window,float> f_factors;
	
	init{
		color<-#white;
		}
	
	action calculate{
		if int(my_room count each.draft)>=1 { 					//L&D flame geometry DRAFT
			flame_tip<-23.9*(wind_speed^(-0.43))*one_of(my_room).mr/(height*width)^0.5-height;
			flame_impigment_range<-0.605*((wind_speed^2/height)^0.22)*(height+width);
		} else {												//L&D flame geometry NO DRAFT
			if width=0 {write "window width is ZERO! "+ self;}
			flame_tip<-12.8*(one_of(my_room).mr/width)-height;
			flame_impigment_range<-height*2/3;
		}
		if balcony {flame_tip<-min([flame_tip,floor_height/2]);} //balconies prevent upwords fires spread
		
		if int(my_room count each.draft)>=1 {
			if downwind {
					flame<-box(flame_impigment_range,width,flame_tip) rotated_by heading;
					flame<-flame translated_by {(flame_impigment_range/2)*cos(heading),(flame_impigment_range/2)*sin(heading),height/3};
					casted_flame<-rectangle(flame_impigment_range,width) rotated_by heading;
					casted_flame<-casted_flame translated_by {(flame_impigment_range/2)*cos(heading),(flame_impigment_range/2)*sin(heading)};
				} else { 
					flame<-nil; //in case of draft conditions, only downwind window radiates
				} 
		} else {
			flame<-box(flame_impigment_range,width,flame_tip) rotated_by heading;
			flame<-flame translated_by {(flame_impigment_range/2)*cos(heading),(flame_impigment_range/2)*sin(heading),height/3};
			casted_flame<-rectangle(flame_impigment_range,width) rotated_by heading;
			casted_flame<-casted_flame translated_by {(flame_impigment_range/2)*cos(heading),(flame_impigment_range/2)*sin(heading)};
		}
		
		if !(first(my_room).draft and upwind ) and flame!=nil {
		loop i over:in_sight{
			float s<-i distance_to self;
			float s_f<-i distance_to self.flame;
			float h_off<-(self.location.z+self.height/2)-(i.location.z+i.height/2); 
			float teta<-self.heading+90-i.heading-90;
			float w_off<-(i distance_to self)*cos(i towards self - self.heading - 90);
			float view_factor<-nil;
			float view_factor_f<-nil;
			//GAS radiation
			if w_off<3/2*self.width{ //if not, configurational factor is zero	
				if w_off<1/2*self.width { //in this case fonfigurational factors are all positive
					
					/*1|2 dividing radiator to 4 sections, see LEE 2009
					 *---
					 *3|4
					 */
					
					float h1<-self.height/2+h_off;
					float w1<-self.width/2+w_off;
					float h2<-h1;
					float w2<-self.width/2-w_off;
					float h3<-self.height/2-h_off;
					float w3<-w1;
					float h4<-h3;
					float w4<-w2;
					
					//calculating view factor
					float a1<-h1/s;
					float b1<-w1/s;
					float a2<-h2/s;
					float b2<-w2/s;
					float a3<-h3/s;
					float b3<-w3/s;
					float a4<-h4/s;
					float b4<-w4/s;
					
					float view_factor1<-(1/2/#pi)*((a1*(1+a1^2)^(-0.5))*#to_rad*atan(b1*(1+a1^2)^(-0.5))+(b1*(1+b1^2)^(-0.5))*#to_rad*atan(a1*(1+b1^2)^(-0.5)));
					float view_factor2<-(1/2/#pi)*((a2*(1+a2^2)^(-0.5))*#to_rad*atan(b2*(1+a2^2)^(-0.5))+(b2*(1+b2^2)^(-0.5))*#to_rad*atan(a2*(1+b2^2)^(-0.5)));
					float view_factor3<-(1/2/#pi)*((a3*(1+a3^2)^(-0.5))*#to_rad*atan(b3*(1+a3^2)^(-0.5))+(b3*(1+b3^2)^(-0.5))*#to_rad*atan(a3*(1+b3^2)^(-0.5)));
					float view_factor4<-(1/2/#pi)*((a4*(1+a4^2)^(-0.5))*#to_rad*atan(b4*(1+a4^2)^(-0.5))+(b4*(1+b4^2)^(-0.5))*#to_rad*atan(a4*(1+b4^2)^(-0.5)));
					view_factor<-view_factor1+view_factor2+view_factor3+view_factor4;
					
				} else { //in this case, two configurational factors are negative
					/*1|2 dividing radiator
					 *---
					 *3|4
					 */
					
					float h1<-self.height/2+h_off;
					float w1<-self.width;
					float h2<-h1;
					float w2<-w_off-self.width/2;
					float h3<-self.height/2-h_off;
					float w3<-w1;
					float h4<-h3;
					float w4<-w2;
					
					//calculating view factor
					float a1<-h1/s;
					float b1<-w1/s;
					float a2<-h2/s;
					float b2<-w2/s;
					float a3<-h3/s;
					float b3<-w3/s;
					float a4<-h4/s;
					float b4<-w4/s;
					
					float view_factor1<-(1/2/#pi)*((a1*(1+a1^2)^(-0.5))*#to_rad*atan(b1*(1+a1^2)^(-0.5))+(b1*(1+b1^2)^(-0.5))*#to_rad*atan(a1*(1+b1^2)^(-0.5)));
					float view_factor2<-(1/2/#pi)*((a2*(1+a2^2)^(-0.5))*#to_rad*atan(b2*(1+a2^2)^(-0.5))+(b2*(1+b2^2)^(-0.5))*#to_rad*atan(a2*(1+b2^2)^(-0.5)));
					float view_factor3<-(1/2/#pi)*((a3*(1+a3^2)^(-0.5))*#to_rad*atan(b3*(1+a3^2)^(-0.5))+(b3*(1+b3^2)^(-0.5))*#to_rad*atan(a3*(1+b3^2)^(-0.5)));
					float view_factor4<-(1/2/#pi)*((a4*(1+a4^2)^(-0.5))*#to_rad*atan(b4*(1+a4^2)^(-0.5))+(b4*(1+b4^2)^(-0.5))*#to_rad*atan(a4*(1+b4^2)^(-0.5)));
					view_factor<-view_factor1+view_factor2+view_factor3+view_factor4;
				}
			}
			
			
			
			//FLAME radiation
			if w_off<3/2*self.width and s_f>0{ //if not, configurational factor is zero	
				if w_off<1/2*self.width { //in this case fonfigurational factors are all positive
					
					/*1|2 dividing radiator
					 *---
					 *3|4
					 */
					
					float h1<-self.height/2+h_off;
					float w1<-self.width/2+w_off;
					float h2<-h1;
					float w2<-self.width/2-w_off;
					float h3<-self.height/2-h_off;
					float w3<-w1;
					float h4<-h3;
					float w4<-w2;
					
					//calculating view factor
					float a1<-h1/s_f;
					float b1<-w1/s_f;
					float a2<-h2/s_f;
					float b2<-w2/s_f;
					float a3<-h3/s_f;
					float b3<-w3/s_f;
					float a4<-h4/s_f;
					float b4<-w4/s_f;
					
					float view_factor1_f<-(1/2/#pi)*((a1*(1+a1^2)^(-0.5))*#to_rad*atan(b1*(1+a1^2)^(-0.5))+(b1*(1+b1^2)^(-0.5))*#to_rad*atan(a1*(1+b1^2)^(-0.5)));
					float view_factor2_f<-(1/2/#pi)*((a2*(1+a2^2)^(-0.5))*#to_rad*atan(b2*(1+a2^2)^(-0.5))+(b2*(1+b2^2)^(-0.5))*#to_rad*atan(a2*(1+b2^2)^(-0.5)));
					float view_factor3_f<-(1/2/#pi)*((a3*(1+a3^2)^(-0.5))*#to_rad*atan(b3*(1+a3^2)^(-0.5))+(b3*(1+b3^2)^(-0.5))*#to_rad*atan(a3*(1+b3^2)^(-0.5)));
					float view_factor4_f<-(1/2/#pi)*((a4*(1+a4^2)^(-0.5))*#to_rad*atan(b4*(1+a4^2)^(-0.5))+(b4*(1+b4^2)^(-0.5))*#to_rad*atan(a4*(1+b4^2)^(-0.5)));
					view_factor_f<-view_factor1_f+view_factor2_f+view_factor3_f+view_factor4_f;
					
					
				} else { //in this case, two configurational factors are negative
					/*1|2 dividing radiator
					 *---
					 *3|4
					 */
					
					float h1<-self.height/2+h_off;
					float w1<-self.width;
					float h2<-h1;
					float w2<-w_off-self.width/2;
					float h3<-self.height/2-h_off;
					float w3<-w1;
					float h4<-h3;
					float w4<-w2;
					
					//calculating view factor
					float a1<-h1/s_f;
					float b1<-w1/s_f;
					float a2<-h2/s_f;
					float b2<-w2/s_f;
					float a3<-h3/s_f;
					float b3<-w3/s_f;
					float a4<-h4/s_f;
					float b4<-w4/s_f;
					
					float view_factor1_f<-(1/2/#pi)*((a1*(1+a1^2)^(-0.5))*#to_rad*atan(b1*(1+a1^2)^(-0.5))+(b1*(1+b1^2)^(-0.5))*#to_rad*atan(a1*(1+b1^2)^(-0.5)));
					float view_factor2_f<-(1/2/#pi)*((a2*(1+a2^2)^(-0.5))*#to_rad*atan(b2*(1+a2^2)^(-0.5))+(b2*(1+b2^2)^(-0.5))*#to_rad*atan(a2*(1+b2^2)^(-0.5)));
					float view_factor3_f<-(1/2/#pi)*((a3*(1+a3^2)^(-0.5))*#to_rad*atan(b3*(1+a3^2)^(-0.5))+(b3*(1+b3^2)^(-0.5))*#to_rad*atan(a3*(1+b3^2)^(-0.5)));
					float view_factor4_f<-(1/2/#pi)*((a4*(1+a4^2)^(-0.5))*#to_rad*atan(b4*(1+a4^2)^(-0.5))+(b4*(1+b4^2)^(-0.5))*#to_rad*atan(a4*(1+b4^2)^(-0.5)));
					view_factor_f<-view_factor1_f+view_factor2_f+view_factor3_f+view_factor4_f;
				}
			}
			self.factors[i]<-view_factor;
			self.f_factors[i]<-view_factor_f;
		}
		}
		
		
		
		calculated<-true;
	}
	
	action impige{
		do change_to_red;
		burn<-true;
		
		if !calculated {
			do calculate;
			list<window> in_my_flame<-(window at_distance 10) where ((each overlaps flame) and (each.location.z<=location.z+height/3+flame_tip) and (each.location.z>=location.z+height/3)) ;
			ask in_my_flame {ask my_room {do ignite;}}
				
			if (casted_flame overlaps trees) and (location.z<=tree_top) {
				burning_area<-burning_area + (casted_flame intersection trees);
			}
		}		
		
		
		
	}
	
	action radiate{
		ask in_sight{
			
			float gas_radiation<-myself.factors[self]*1*((one_of(myself.my_room).Tr)^4-normal_temperature^4)*(5.76*10^(-8));
			float flame_radiation<-myself.f_factors[self]*(1-exp(-0.3*myself.flame_impigment_range))*(813^4-normal_temperature^4)*(5.76*10^(-8));
			do recive_radiation(gas_radiation+flame_radiation);
		}
	}
	
	action recive_radiation(float influx){
		switch influx{
			match_between [0,12500] {radiation_time<-0;}
			match_between [12500,15000] {radiation_time<-radiation_time + step;}
			match_between [15000,17500] {radiation_time<-radiation_time + step*30/25;}
			match_between [17500,20000] {radiation_time<-radiation_time + step*30/10;}
			match_between [20000,30000] {radiation_time<-radiation_time + step*30/7;}
			match_between [30000,3000000] {radiation_time<-radiation_time + step*30/1;}
		}
		if radiation_time>30*60 {
			ask first(my_room) {do ignite;}
			ignited_by_rad<-true;
			time_to_rad<-time;
		}
	}
	
	action recive_rad_from_veg (float timing){
		if !my_room[0].burn and !my_room[0].burnt_down {ignition_time_by_veg_rad<-time+timing;}
	}
	
	aspect over_view{
		draw box(0.5,width,height) transformed_by {heading,1} color: color ;
		if balcony {draw box(1,width,0.1) transformed_by {heading,1} translated_by{0,0,floor_height/2} color:#black;}
	}
	aspect base_3D{
		if upwind {color<-#green;}
		if downwind {color<-#yellow;}
		draw box(0.5,width,height) transformed_by {heading,1} color: color ;
		draw flame color:#gold;
		if balcony {draw box(1,width,0.1) transformed_by {heading,1} translated_by{0,0,floor_height/2} color:#black;}
		if show_in_sight{
			loop i over:in_sight{
				draw line([self,i]) color:#blue;  //shows line of sight to other windows
			}
		}
	}
	aspect base{
		draw rectangle(0.5,width) transformed_by {heading,1} color: color border:#grey;
		draw casted_flame color:#gold;
		}
}

species door parent:element{
	list<room> my_rooms<-nil;
	
	init{
		color<-#blue;			
	}
	
	action impige{
		if int(my_rooms count each.burn)>=1 {ask (my_rooms - (my_rooms where each.burn)) {do ignite;}}
	}
	
	aspect over_view{
		draw cube(0.8)  color: color depth:floor_height;
	}
	
	aspect base{
		draw rectangle(0.8,0.8) color: color;
		}
}

species tree parent:element{
	init{
		color<-#green;	
	}
	aspect base{
		draw shape color: color;
		}
}

species brand {
	rgb color;
	init{
		color<-#red;
	}
	
	action ignite{
		if location overlaps world {
			geometry my_projection <- box(0.1,0.1,100) at_location self.location;
			list<room> landed_on<-(room where(each overlaps my_projection)) where each.flame_roof;

			if length(landed_on)>0 {ask landed_on {
				do ignite;
				}
			}
			if flip(brand_veg_ignition_factor) and self overlaps trees{
				burning_area<-burning_area + (geometry(self)+1);
			}
		}
		
	}
	
	aspect base{
		draw shape color: color;
		}
}


experiment serial_apts type:batch  repeat:1 until:time>=time_limit {
	parameter var:time_limit <-1*#hour;
	parameter "serial apt" category:"simulation" var:go_serial<-true;
	parameter category:"simulation" var:ignitions<-1; 
	parameter "step duartion (secondes)" category:"simulation" var:step<-60;
	parameter "create trees?" category:"simulation" var:create_trees<-true;
	parameter "floor_height" var:floor_height<-2.8;
	parameter "balcony ratio" category:"window" var:balcony_ratio<-0.087;
	parameter "trees buffer" category:"trees" var:trees_buffer<-0;
	parameter "trees burn range" category:"trees" var: trees_burn_range<-75 #m;
	parameter "wind direction" category:"meteorology" var:wind<-90; //among:[45,135,180,225,270,315];
	parameter "wind speed" category:"meteorology" var:wind_speed among: [4.0,6.0,8.0];
	parameter "vegetation type" category:"trees" var:veg_type<-"shrub";// among:["shrub","trees_5m","trees_10m"];
	parameter "vegetation fron burning time (s)" category:"trees" var:veg_front_burn_time<-120#s ;//among:[300#s,240#s,180#s,120#s];
	parameter "barnding probability" category: "fire" var:brand_prob<-0.1;
	parameter "impigement range" category: "fire" var:impigement_range<-0.5;
	parameter "fire load (wood kg/m2)" category:"fire" var:L_general<-16;
	parameter "brand veg ignition factor" category:"fire" var:brand_veg_ignition_factor<-1/4;
	parameter "output file" category:"output" var:output_file_name<-"../doc/results.csv";
	parameter var:polygons_recording<-false;
	parameter "CRS" category:"GIS" var:CRS<-"EPSG:2039"; //do not use unprojected coordinates
	parameter "buildings footprints" category:"GIS" var:shape_buildings_file<-shape_file("../includes/city/city_buildings1.shp",CRS);
	parameter "trees" category:"GIS" var:shape_trees_file<-shape_file("../includes/city/city_trees_fine.shp",CRS);
	parameter "bounds" category:"GIS" var:shape_bounds_file<-shape_file("../includes/city/city_bounds.shp",CRS);
	parameter "rooms" category:"GIS" var:shape_rooms_file<-shape_file("../includes/new_city/rooms1.shp",CRS);
	parameter "doors" category:"GIS" var:shape_doors_file<-shape_file("../includes/new_city/doors.shp",CRS);
	parameter "windows" category:"GIS" var:shape_windows_file<-shape_file("../includes/new_city/windows.shp",CRS);
	parameter category:"simulation" var:report_init<-false;
	//serial run:
	parameter "to_explore" var:to_explore min:0 max:1613 step:1; //max should be the (number of ground floor apartments -1) or (number of buildings -1)
	
	action _step_ {
		ask simulations{
			write ""+self+" run took: " + ((machine_time-zero_time)/1000);
			//write output:
			write ["MME model - ","simulated time:",(time),"PARAMETERS:",step,tune_factor,rooms_file,shape_trees_file.name,balcony_ratio,L_general,ignitions,to_explore,origin_of_fire,origin_of_fire_building,origin_of_fire_floor,wind_speed,wind,burning_speed,create_trees,trees_burn_range,impigement_range,tree_top,trees_buffer,veg_type,veg_front_burn_time,brand_prob,brand_veg_ignition_factor,"RESULTS :",veg_was_ignited,time_to_veg,brand_were_released,time_to_brands,ignited_by_rad,time_to_rad,time_to_first_FO,time_to_next_apt,time_to_next_floor,next_ignited_building,time_to_next_building,d_to_next_building,angel_to_next_building,"MAXS :",max(burnt_area_series),max(burning_room_series),max(burning_apt_series),max(burning_floor_series),top_burning_floor,max(burning_buildings_series)];
			//save output:
			save(["MME model - ","simulated time:",(time),"PARAMETERS:",step,tune_factor,rooms_file,shape_trees_file.name,balcony_ratio,L_general,ignitions,to_explore,origin_of_fire,origin_of_fire_building,origin_of_fire_floor,wind_speed,wind,burning_speed,create_trees,trees_burn_range,impigement_range,tree_top,trees_buffer,veg_type,veg_front_burn_time,brand_prob,brand_veg_ignition_factor,"RESULTS :",veg_was_ignited,time_to_veg,brand_were_released,time_to_brands,ignited_by_rad,time_to_rad,time_to_first_FO,time_to_next_apt,time_to_next_floor,next_ignited_building,time_to_next_building,d_to_next_building,angel_to_next_building,"MAXS :",max(burnt_area_series),max(burning_room_series),max(burning_apt_series),max(burning_floor_series),top_burning_floor,max(burning_buildings_series)]) type: "csv" to: output_file_name rewrite:false;
			/*allows saving fire spread in vegetation*/ if veg_was_ignited and polygons_recording {save projected_burn_area type:"shp" to:"../doc/city/tuned_proj/projections_building_case_"+to_explore+".shp"  with:[timing::'timing'];}
		}	
	}
	
	
}


	
experiment regular_display type:gui {
	parameter var:time_limit <-1*#hour;
	parameter category:"simulation" var:report_init<-false;
	parameter category:"simulation" var:ignitions<-1; 
	parameter "serial" category:"simulation" var:serial<-0;
	parameter "step duartion (secondes)" category:"simulation" var:step<-60;
	parameter "create trees?" category:"simulation" var:create_trees<-false; //change to 'true'
	parameter "output file" category:"simulation" var:output_file_name<-"../doc/results.csv";
	parameter "floor_height" var:floor_height<-2.8;
	parameter "balcony ratio" category:"window" var:balcony_ratio<-0.087;
	parameter "trees buffer" category:"trees" var:trees_buffer<-0;
	parameter "vegetation type" category:"trees" var:veg_type<-"shrub" among:["shrub","shrub-dry","trees_5m","trees_10m"];
	parameter "vegetation fron burning time (s)" category:"trees" var:veg_front_burn_time<-120#s among:[300#s,240#s,180#s,120#s];
	parameter "trees burn range" category:"trees" var: trees_burn_range<-75 #m;
	parameter "wind direction" category:"meteorology" var:wind<-90;
	parameter "wind speed" category:"meteorology" var:wind_speed among:[4.0,6,8];
	parameter "barnding probability" category: "fire" var:brand_prob<-0.1;
	parameter "impigement range" category: "fire" var:impigement_range<-0.5;
	parameter "fire load (wood kg/m2)" category:"fire" var:L_general<-16;
	parameter "brand veg ignition factor" category:"fire" var:brand_veg_ignition_factor<-1/4;
	parameter "CRS" category:"GIS" var:CRS<-"EPSG:2039"; //do not use unprojected coordinates
	parameter "buildings footprints" category:"GIS" var:shape_buildings_file<-shape_file("../includes/algorithm/sample_building_lines.shp",CRS);
	// not trees file for the algorithm sample. Use your own: parameter "trees" category:"GIS" var:shape_trees_file<-shape_file("../includes/city/city_trees_fine.shp",CRS);
	parameter "bounds" category:"GIS" var:shape_bounds_file<-shape_file("../includes/sample_bounds.shp",CRS);
	parameter "rooms" category:"GIS" var:shape_rooms_file<-shape_file("../includes/algorithm/sample_rooms.shp",CRS);
	parameter "doors" category:"GIS" var:shape_doors_file<-shape_file("../includes/algorithm/sample_doors.shp",CRS);
	parameter "windows" category:"GIS" var:shape_windows_file<-shape_file("../includes/algorithm/sample_windows.shp",CRS);
	output{
		
		display fast_over_view type:opengl {
			species building_line aspect:base_3D;
			
			species window aspect: base_3D;
			species roof_flame aspect: base_3D;
			graphics name:"burnin trees"{
			
				draw burning_area  color:#pink ;
				
				
				draw trees depth: tree_top empty:true border: #green color:#green;	
			}
		}
		
		
		
		
		display development {
			chart "burning/burnt down" type: series size: {1,1/2} position: {0, 0} legend_font_size:15 label_font_size:15 tick_font_size:15{
				 datalist ["apts","floors","buildings"] value:[last(burning_apt_series),last(burning_floor_series),last(burning_buildings_series)] color:[#blue,#green,#brown];
			}
			chart "rooms" type: series size: {1,1/2} position: {0, 1/2} legend_font_size:15 label_font_size:15 tick_font_size:15{
				 datalist ["rooms","burnt down rooms"] value:[last(burning_room_series),length(room where (!(each.burn) and (each.burnt_down)))] color:[#red,#black];
			}
		} 
		
		
		display map  {	
			graphics name:"ranges"{
				draw impingemnt_area color:#purple empty:true;
				draw trees empty:true border: #green color:#green;
			}
			agents name:"rooms" value:room where (each.floor=1) aspect: base ;
			species window aspect: base ;
			agents name:"doors" value:door where (each.floor=1) aspect: base  ;
			species brand aspect:base;
			species projected_burn_area aspect:base;
			graphics name:"burnin trees" transparency: 0.5{
			
				
				draw veg_radiation_area_60s color:#red;
				draw veg_radiation_area_120s color:#darkred;
				draw veg_radiation_area_180s color:#gold ;
				draw veg_radiation_area_300s color:#brown;
				draw burning_area color:#black;
				
				
			}
		}
		
	}
}

