/**
 *  Fire Spread Model for Mediterranean and Middle Eastern Cities - Data Generation Algorithm
 *  Authors: Yonatan Shaham (yjonas83@gmail.com), Itzhak Benenson (bennya@post.tau.ac.il)
 *  Licence: GPLv3.0 see: http://www.gnu.org/licenses/gpl-3.0.txt 
 * 
 * This implimentation does not include special cases.
 */

model MME_alg

global torus:false {
	//global
	float step <- 1 #mn;
	float sim_start;
	float L_general<-16; //mean fuel load per sq.m.
	
	//GIS
	string buildings_file<-"../includes/sample_buildings.shp";
	string CRS<-"EPSG:2039";
	file shape_buildings_file<-shape_file(buildings_file,CRS);
	file shape_bounds_file<-shape_file(bounds_file,CRS);
	
	string bounds_file<-"../includes/sample_bounds.shp";
	string file_prefix<-"sample_";
	geometry shape <- envelope(shape_bounds_file);	
	
	bool create_apts<-true;
	float Aa<-80; // typical apartment area in sq.m.
	float Ra<-12.25; // typical room area in sq.m.
	float Rw<-Ra^0.5;
	float floor_height<-2.8;
	
	float w_h<-1; //defualt window height in m.
	float w_w<-1; //defualt window width in m.
	float W<-2; //Minimal length of a wall with a window
	int outer_rooms_in_apt<-3;
	
	bool reoprt_times<-true;
	float start_time;
	float step_time;
	
	
	building_line b;
	int b_number<-0;
	geometry outside;
	
	geometry temp0;
	geometry temp1;
	geometry best_mbr;
	geometry best_line;
	point ancor;
	point farthest;
	point second;
	
	int number_of_sections;
	list<room> outer_ring;
	list<room> current_ring;
	list<room> current_corridor;
	
	init{
		start_time<-machine_time;
		create building_line from:shape_buildings_file with:[BLDG_HT::float(read("BLDG_HT")),L::float(read("L")),flame_roof::int(read("roof"))];
		outside<-world.shape;
		loop i over:building_line{
			i.UNIQ_ID<-i;
			i.floors<-max([1,floor(i.BLDG_HT/floor_height)+1]);
			outside<-outside-i.shape;
		}
		
		
		write "Starting MME data generation algorithm:";
		write "One building will be handeled in each simulation step";
		write "Corridors are colored brown.";
		write "=======================================================";
	}
	
	reflex go{
		if b_number>length(building_line)-1 { write "DONE"; do pause;}
		step_time<-machine_time;
		
		ask temp_wall {do die;}
		ask poly {do die;}
		ask temp_room {do die;}
		ask wall {do die;}
		ask control_dot {do die;}
		
		b<-building_line[b_number];
		int apt<-1;
		
		write "now working on "+b;	
		
		float foot_print_area<-b.shape.area;
		
		list<point> my_points<-b.shape.points;
				
		geometry con<-convex_hull(b.shape);
		
		list<geometry> con_lines<-nil;
		loop i from:0 to: length(con.points)-2{
					con_lines<-con_lines + line([con.points[i],con.points[i+1]]);
					}
		float mbr_area<-100000;
		best_mbr<-nil;
		
		geometry max_line<-nil;
		loop i over:con_lines{
			float teta<-(first(i.points) towards last(i.points));
			geometry temp_con<-con rotated_by (-teta);
			geometry temp_mbr<-envelope(temp_con);
			
			if temp_mbr.area<mbr_area {
				mbr_area<-temp_mbr.area;
				
				best_mbr<-(temp_mbr) rotated_by teta;
				point move<-{(b.shape.points with_max_of each.x).x - (best_mbr.points with_max_of each.x).x,(b.shape.points with_max_of each.x).y - (best_mbr.points with_max_of each.x).y};
				best_mbr<-(best_mbr translated_by move);
				best_line<-i;
				max_line<-i;
			}
			
		}
		temp0<-best_mbr;
		list<geometry> mbr_lines<-nil;
		loop i from:0 to: length(best_mbr.points)-2{
					mbr_lines<-mbr_lines + line([best_mbr.points[i],best_mbr.points[i+1]]);
					}
		list<geometry> long_sides<-mbr_lines where (each.perimeter = (mbr_lines with_max_of each.perimeter).perimeter);
	
		
		
		best_line<-long_sides with_min_of (each.perimeter - (mbr_lines closest_to each).perimeter );
		int p<-0;
		
		ancor<-best_line.points with_max_of each.y;
		second<-best_line.points with_min_of each.y;
		farthest<-(best_mbr.points) with_max_of (each distance_to ancor);
		int move1_times<-((ancor distance_to second)/Rw)+1;
		int move2_times<-min([((ancor distance_to farthest)/Rw)+1,((second distance_to farthest)/Rw)+1]);
		float heading<-ancor towards second;
				
		point out_vector<-nil;
		float out_size<-( (Rw/2)^2 + (Rw/2)^2 )^0.5;
		
		out_vector<-{cos(heading-45)*(out_size),sin(heading-45)*(out_size)};
		point move_vector<-{cos(heading)*(Rw),sin(heading)*(Rw)};
		point move_vector2<-{cos(heading+90)*(Rw)+0.01,sin(heading+90)*(Rw)+0.01};
		if farthest.x < second.x {
			move_vector2<-{cos(heading-90)*(Rw)+0.01,sin(heading-90)*(Rw)+0.01};
			out_vector<-{cos(heading+45)*(out_size),sin(heading+45)*(out_size)};
		}
		p<-0;
		int t<-1;
		point normal_point<-(best_line closest_points_with farthest)[0];
		geometry normal<-line([farthest,normal_point]);
		
		
		
		loop times: move2_times {
			p<-0;
			loop times: move1_times{
				create poly number:1{
					shape<-rectangle(Rw,Rw) rotated_by heading;
					location<-ancor;
					location<-location+out_vector;
					location<-location+move_vector*p;
					location<-location+move_vector2*t;
					UNIQ_ID<-b.UNIQ_ID;
					section<-floor(p/5);
					section_x<-1+p-section*5;
					section_y<-t;
					floors<-b.floors;
				}
				p<-p+1;
			}
			t<-t+1;
		}
			
		loop i over:poly{
			building_line my_building <- one_of((building_line at_distance 10 where (each.UNIQ_ID=i.UNIQ_ID)) overlapping i);
			geometry temp<-i.shape intersection my_building.shape;
			if !(temp=nil){
				create temp_room number:1 {
					UNIQ_ID<-my_building.UNIQ_ID;
					shape<-temp;
					section<-i.section;			
					section_x<-i.section_x;
					section_y<-i.section_y;
					floors<-b.floors;
				}
			}
		}
		
		ask (temp_room sort_by each.shape.area){
			if shape.area<7{
				if shape.area<1.5 {do die;} //rooms of this size area usually very small slivers
				temp_room best<-nil;
				
				loop i from:0 to: length(shape.points)-2{
					create temp_wall number:1 {
						shape<-line([myself.shape.points[i],myself.shape.points[i+1]]);
					}
				}
				
				ask temp_wall{
					my_rooms<-(temp_room at_distance 10) overlapping (circle(0.5) at_location self.location);	
					if length(my_rooms)<2 {do die;}
				}
				
				temp_wall best_wall<-temp_wall with_max_of each.shape.perimeter;
				best<- one_of(best_wall.my_rooms-self);
				if best!=nil{
					ask best{
						shape<-shape union myself.shape;
						shape<-shape+0.1;
						shape<-shape-0.1;
						
						section_x<-min([myself.section_x,best.section_x]);
						section_y<-min([myself.section_y,best.section_y]);
					}	
				} 
				
				do die;
			}
		}
		
		ask temp_room{
			create room number:1 {	
				floors<-myself.floors;
				UNIQ_ID<-myself.UNIQ_ID;
				shape<-myself.shape.contour;
				shape<-polygon(remove_duplicates(shape.points));
				section<-myself.section;			
				section_x<-myself.section_x;
				section_y<-myself.section_y;
				
				number_of_sections<-max([number_of_sections,section]);	
			}		
		}
		write "number_of_sections: "+number_of_sections;
		//walls creation
		ask room where (!each.done){
			loop i from:0 to: length(shape.points)-2{
				create wall number:1 {
					shape<-line([myself.shape.points[i],myself.shape.points[i+1]]);
					floor<-1;
					heading<-shape towards myself.shape;
					myself.my_walls<-myself.my_walls+self;
				}
			}
		}
		
		
		ask room where (!each.done){
			if (self.shape+0.2) overlaps outside{
				outer_ring<-outer_ring+self;
			}
			
			list<room> near<-(room at_distance 10) overlapping (self.shape+1);
			
			near<-remove_duplicates(near);
			loop i over:near{
				if i.section_x=section_x or i.section_y=section_y{
					create control_line number:1{
						shape<-line([myself.location,i.location]);
					}
				}
			}
		}
		
		ask wall where (!each.done){
			list<control_line> inters<-control_line overlapping self;
			loop i over:inters{
				create control_dot{
					location<-i inter myself;
				}
			}
		}
		
		ask control_dot{
			
			near<-(control_dot at_distance 3) overlapping (self+1);
			near<-near-self;
		}
		
		
		loop i over:control_dot{
			
			ask i.near {
				
				do die;
			}
		}
		
		//MERGING SECTIONS WITH MAX_X<3
		write "MERGING SECTIONS WITH MAX_X<3";
		loop section from:0 to: number_of_sections{
			write "===============================";
			write "working on section: "+section;
			list<room> in_section<-room where(!each.done and each.section=section);
	
			int max_y<-(in_section with_max_of each.section_y).section_y;
			int max_x<-(in_section with_max_of each.section_x).section_x;
	
			if max_x<3 and section>0 and section=number_of_sections{
				list<room> in_prev_section<-room where(each.UNIQ_ID=b.UNIQ_ID and each.section=section-1);
				
				int prev_max_y<-(in_prev_section with_max_of each.section_y).section_y;
				int prev_max_x<-(in_prev_section with_max_of each.section_x).section_x;
		
				ask in_section{

					section<-section-1;
					section_x<-section_x+prev_max_x;
	
				}
				
			}
		}
		
		
		write "END - MERGING SECTIONS WITH MAX_X<3";
		write "===============================";
		write "START - SECTION WORK";
		loop section from:0 to: number_of_sections{
			write "===============================";
			write "working on section: "+section;
			list<room> in_section<-room where(!each.done and each.section=section);
			write "in_section: "+in_section;
			int max_y<-(in_section with_max_of each.section_y).section_y;
			int max_x<-(in_section with_max_of each.section_x).section_x;
			write "max x: "+max_x+" max y: "+max_y;
			
			if max_x>=3 and max_x<=7 {
				int corridor_x<-int(max_x/2)+1;
				if max_x=4 {corridor_x<-2;}
		
				list<room> to_corridor<-in_section where (each.section_x=corridor_x and each.section_y<length(in_section where (each.section_x=corridor_x)));
				write to_corridor accumulate (each.name + " section_x: "+each.section_x+" section_y: "+each.section_y);
				list<geometry> corridor_shape<-to_corridor accumulate each.shape;
				point corridor_center;
				create room number:1 {
					shape<-union(corridor_shape);
					shape<-shape+0.1;
					shape<-shape-0.1;
					corridor<-true;
					corridor_center<-self;
					
				}
				ask to_corridor{
					
					do die;
				}
				in_section<-room where(!each.done and each.section=section);
				
				//OUTERRING
				
				geometry section_shape<-union(in_section accumulate (each.shape+0.1));
				outside<-world.shape-section_shape;
				outside<-outside+1;
				current_ring<-nil;
				
				current_ring<-in_section where (!each.corridor and (each overlaps outside)) ;
				write current_ring;
				
				
				
				write corridor_x;
				ask current_ring{
				
				}
				room current_room<-first(current_ring where (each.section_x=(corridor_x-1) and each.section_y=1));
				
				list<room> sorted_ring<-nil;
				sorted_ring<-sorted_ring+current_room;
				current_ring<-current_ring-current_room;
				room next_room<-nil;
				list<room> pot;
				
				bool continue<-true;
				loop while: length(current_ring)>0 and continue{
					pot<-current_ring where (each overlaps (current_room+0.5));
					if length(pot)=0 {continue<-false;}
					next_room<-pot with_max_of (each distance_to corridor_center);
					
					sorted_ring<-sorted_ring+next_room;
					current_ring<-current_ring-next_room;
					current_ring<-remove_duplicates(current_ring);
					
					current_room<-next_room;
				}
				
				current_ring<-sorted_ring;
				int outer_in_current_apt<-0;
				loop i over: current_ring{
					
					if outer_in_current_apt<outer_rooms_in_apt{
						i.apt<-apt;
						outer_in_current_apt<-outer_in_current_apt+1;
					} else {
						apt<-apt+1;
						i.apt<-apt;
						outer_in_current_apt<-1;
					}
					
				}
				ask current_ring {color<-#white;}
				//INNER RING
				
				current_ring<-in_section where (!each.corridor and !(each overlaps outside)) ;
			
				ask current_ring {
					color<-#pink;
					list<room> next_to_me<-(room where (each.apt!=0 and each.section=self.section)) overlapping (self+0.5);
					
					int common_apt<-0;
					int max_counter<-0;
					loop j from:0 to:max(next_to_me accumulate each.apt){
						int counter<-length(next_to_me where (each.apt=j));
						
						if counter>max_counter{
							max_counter<-counter;
							common_apt<-j;
							
						}
					}
					
					self.apt<-common_apt;
				}
			} else {
				//special cases
				ask in_section{apt<-100;}
			}
			
			
			
			
			if length(in_section)>0{
				
				
				list<int> apts_in_section<-in_section accumulate each.apt;
				list<room> to_re_apt;
				loop j from:0 to:max(apts_in_section){
					int rooms_in_apt<-length(apts_in_section where (each=j));
					if rooms_in_apt<3 {
						ask in_section where (each.apt=j and !each.corridor){
							self.apt<-0;
							color<-#red;
							to_re_apt<-to_re_apt+self;
						}
					}
				}
				
				ask to_re_apt{
					color<-#pink;
					list<room> next_to_me<-(room where (each.apt!=0 and each.section=self.section)) overlapping (self+0.5);
					
					int common_apt<-0;
					int max_counter<-0;
					loop j from:0 to:max(next_to_me accumulate each.apt){
					int counter<-length(next_to_me where (each.apt=j));
					
					if counter>max_counter{
						max_counter<-counter;
						common_apt<-j;
						
						}
					}
					
					self.apt<-common_apt;
				}
				
				apts_in_section<-remove_duplicates(in_section accumulate each.apt);
				loop i from:min(apts_in_section) to: max(apts_in_section){
					list<room> in_apt<-in_section where (each.apt=i);
					list<room> in_apt_with_corridor;
					ask in_apt{
						list<room> my_rooms<-in_section where (each.corridor and (each overlaps (self+0.1)));
						if length(my_rooms)>0 {in_apt_with_corridor<-in_apt_with_corridor+self;}
					}
					ask in_apt_with_corridor closest_to best_line {color<-#gold;}
				}
				
				
				ask control_dot{
					my_rooms<-(room at_distance 10) overlapping (self);
					
					my_rooms<-my_rooms where (each.apt>0);
					my_rooms<-my_rooms where (each.section=section);
					if length(my_rooms)=2{
						if my_rooms[0].apt=my_rooms[1].apt {
							create door number:1{
								location<-myself;
							}
						}
					}
				}
				
				ask control_dot{
					my_rooms<-(room at_distance 10) overlapping (self);
					
					my_rooms<-my_rooms where (each.apt=0 or each.color=#gold);
					my_rooms<-my_rooms where (each.section=section);
					if length(my_rooms)=2{
						if true {
							create door number:1{
								location<-myself;
							}
						}
					}
				}
				
				
				
				section<-section+1;
			}
		
		}
		
		
		ask room where (!each.done) {
			if !corridor {
				if b.L!=0 {Ltag<-b.L;} else {Ltag<-L_general;}
			} else {
				Ltag<-0; color<-#brown;
			}
		}
		//WINDOWS
		outside<-world.shape-b.shape;
		outside<-outside+Rw*0.35;
		ask wall where (each.location overlaps outside and each.shape.perimeter>=W  ) {
			
				if shape.perimeter<2*W {
					loc<-self;
					create window {
						location<-myself.loc; floor<-myself.floor;  height<-myself.my_window_height; width<-myself.my_window_width;Ad<-height*width; balcony<-myself.balcony;
						heading<-myself.shape.points[0] towards myself.shape.points[1]+90;
					}
				} else {
					list<point> window_points;
					float teta<-shape.points[0] towards shape.points[1];
					loop i from:1 to: int(shape.perimeter/W){
						window_points<-shape.points[0]+{cos(teta)*W*i,sin(teta)*W*i};
					}
					loop i over:window_points {
						create window {location<-i; floor<-myself.floor;  height<-myself.my_window_height; width<-myself.my_window_width;Ad<-height*width; balcony<-myself.balcony;}
					} 
				}
				
		}
		ask window where (!each.done) {
			my_room<-(room where (!each.done) ) overlapping (self+0.5);
			if length(my_room)>0{
				if my_room[0].corridor {
					write ""+self+" of corridor";
					do die;
				}	
				if length(my_room)>1{
					write ""+self+" inside building";
					do die;
				}
			}
		}
		ask window where (!each.done) {
			my_room<-(room where (!each.done) ) overlapping (self+0.5);
			if ! ((location+{cos(heading)*2,sin(heading)*2}) overlaps outside){
				heading<-heading+180;
			}
			done<-true;
		}
		
		ask door where (!each.done) {id<-door index_of self;}
		ask room where (!each.done) {
			if shape.area<1.5 {do die;} //rooms of this size area usually very small slivers
		}
		ask wall {done<-true;}
		ask room {done<-true;}
		ask b {done<-true;}
		ask door {done<-true;}
		ask window {done<-true;}
		save room type:"shp" to:"../includes/algorithm/"+file_prefix+"rooms.shp" with:[floors::'floors',apt::'apt',UNIQ_ID::"UNIQ_ID",corridor::"corridor",Ltag::"L"];
		save window type:"shp" to:"../includes/algorithm/"+file_prefix+"windows.shp" with:[width::"width",height::"height",heading::"heading"];
		save door type:"shp" to:"../includes/algorithm/"+file_prefix+"doors.shp" with:[id::'id'];
		save building_line  where (each.done) type:"shp" to:"../includes/algorithm/"+file_prefix+"building_lines.shp" with:[floors::'floors',UNIQ_ID::"UNIQ_ID",flame_roof::"roof"];
		
		b_number<-b_number+1;
		if reoprt_times {write "this step took :" + ((machine_time-step_time)/1000); write "total runtime: " + ((machine_time-start_time)/1000);}
		write "---";
	}

}

species control_line{
	aspect base{
		draw self color:#blue;
	}
}

species temp_wall{
	list<temp_room> my_rooms;
}

species poly{
	int apt;
	int UNIQ_ID;
	int section;
	int section_x;
	int section_y;
	int floors;
	aspect base{
		draw shape color:#purple empty:true;
		draw ""+section+",x:"+section_x+"y:"+section_y font: font('Default', 11, #bold) color:#black at:location-{2,0,0};
	}
	aspect paper{
		draw shape color:#purple empty:true border:#purple;
	}
}

species building_line{
	int UNIQ_ID;
	float BLDG_HT;
	float L;
	rgb color<-#white;
	bool done;
	int floors;
	int flame_roof;
	aspect base{
		draw shape color:color border:true;
	}
	
}

species temp_room{
	int apt;
	int UNIQ_ID;
	int floors;
	int section;
	int section_x;
	int section_y;
	
	aspect base{
		draw shape color:#gold empty:true;
		draw ""+section+",x:"+section_x+"y:"+section_y font: font('Default', 11, #bold) color:#black at:location-{2,0,0};
	}
}

species control_dot {
	rgb color<-#blue;
	int floors;
	list<room> my_rooms;
	int UNIQ_ID;
	bool ignited;
	list<control_dot> near;
	
	init{
		shape<-circle(0.1);
	}
	
	aspect base {
		draw shape color: color empty:true;
	} 
}

species element{
	rgb color <- #white;
	int floor<-nil;
	
	bool done;
	
	aspect base {
		draw shape color: color;
				
	}
}

species room parent:element{
	
	list<window> my_windows<-nil;
	list<door> my_doors<-nil;
	list<wall> my_walls<-nil;
	int UNIQ_ID;
	list<wall> temp_walls;
	float BLDG_HT<-nil;
	int apt;
	int floors;
	bool flame_roof<-false; 	//is the roof flamable ?
	float Ltag<-L_general #kg/#m2; 		//fuel load density
	int section;
	int section_x;
	int section_y;
	bool corridor;
	bool done<-false;
	bool inside;
	string apt_UNIQ;
	bool outer<-false;
	int ring_index;
	
	aspect base {
		draw shape color: color border:#black;
		draw ""+section+",x:"+section_x+"y:"+section_y font: font('Default', 11, #bold) color:#black at:location-{2,0,0};
	}
	
	aspect ring{
		draw shape color: color border:#black;
		draw ""+ring_index font: font('Default', 14, #bold) color:#black at:location;
	}
	
	aspect apt{
		draw shape color: color border:#black;
		draw ""+apt font: font('Default', 14, #bold) color:#black at:location;
	}
	aspect paper{
		draw shape color: #white border:#black;
		draw ""+apt font: font('Default', 14, #bold) color:#black at:location;
	}
}

species wall parent:element{ //walls are used only for creating doors and windows
	int i_window<-0;
	int i_door<-0;
	int is_door<-nil;
	float heading<-nil;
	
	float my_window_height<-w_h;
	float my_window_width<-w_w;
	bool balcony<-false;
	room my_room;
	
	//variable for room creation algorithm
	list<room> temp_rooms;
	point loc;
	
	action create_windows_and_doors{
		if !(my_room=nil){
			
		}
		
		if i_window=1{
				
			}
		if i_door=1{
			
			}
	}
	
	aspect base{
		draw (shape+0.1) color: #black empty:true;
	}
}

species window parent:element{
	float heading;
	list<room> my_room<-nil;
	
	float internal_heading<-nil;
	float heading_check<-nil;
	float width<-nil;
	float height<-nil;
	float Ad<-nil;
	bool balcony<-false;
	int floors;
	
	init{
		color<-#white;
		}
	
	aspect over_view{
		draw box(0.5,width,height) transformed_by {heading,1} color: color ;
	}
	
	aspect base{
		draw rectangle(0.5,width) transformed_by {heading,1} color: color border:#black;
	}
}

species door parent:element{
	list<room> my_rooms<-nil;
	int floors;
	int id;
	bool done;
	init{
		color<-#blue;			
	}
		
	aspect base{
		draw rectangle(0.8,0.8) color: color;
	}
}
	
experiment running_the_algorithm type:gui {
	
	parameter var:CRS<-"EPSG:2039"; //do not use unprojected coordinates
	parameter "buildings footprints" var:shape_buildings_file<-shape_file("../includes/sample_buildings.shp",CRS);
	parameter "bounds" var:shape_bounds_file<-shape_file("../includes/sample_bounds.shp",CRS);
	parameter "file prefix" var:file_prefix<-"sample_";
	parameter "create apartments?"  var:create_apts<-true;
	parameter "typical apartment area (sq.m.)" var: Aa<-90; 
	parameter "typical room area (sq.m.)" var: Ra<-12.25; 
	parameter "defualt window height (m)" var: w_h<-1; 
	parameter "defualt window width (m)" var:  w_w<-1; 
	parameter "floor height (m)" var:floor_height<-2.8;
	parameter "fire load (wood kg/sq.m.)" var:L_general<-16;
	
	output{
		display apt{
			graphics "current building mbr" {
				draw outside color:#grey ;
			}			
			species building_line aspect:base;
			species room aspect: apt;
			species window aspect:base;
			species door aspect:base;
		}
	}
}

