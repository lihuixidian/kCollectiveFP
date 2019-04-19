// PLS.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include <chrono>//计时
#include <iostream>
#include <fstream>
#include <sstream>

#include "checkin_dump.h"
#include "query_generator.h"
#include "probability_functions.h"
#include "venue_checkin_dump.h"
#include "pino.h"



#define RELEASE_VER 1

#define CANDIDATES_GENERATION 1     //candidate generates
#define GEN_FROM_10162 1       // From 10162 (1) or 2321 (0) dataset
#define PICK_FROM_UNIQUE 1   // In CANDIDATES_GENERATION, pick from unique coordinates (set) but not all check-in logs (vector).
#define CHECKINS_EXCLUDING 1  // Exclude check-ins that appear in candidates.
#define DATA_LOADING 1        




int _tmain(int argc, _TCHAR* argv[])
{

	auto begin_time = std::chrono::high_resolution_clock::now();
	auto end_time = std::chrono::high_resolution_clock::now();
	
	checkin_dump data_dump;
	data_dump.set_string_checkin_meta_format(true, false, checkin::citf_hour_min, '\t', ',', '-', ' ', ':', 5, 0, 2, 0, 0, 3);

	venue_checkin_dump vc_data_dump;

	std::vector<geo_coordinate> candidate_vector;
	unsigned candidate_count = 600; //the number of total candidates

		

#if CANDIDATES_GENERATION
		std::cout << ">>>>> CANDIDATES_GENERATION >>>>>\n" << std::endl;

		std::vector<geo_coordinate> coordinate_vector;

		std::cout << "1  Loading check-in data (no time) from checkins.txt to a geo_coordinate vector." << std::endl;
		begin_time = std::chrono::high_resolution_clock::now();
#if GEN_FROM_10162
		data_dump.load_string_checkin_without_time_to_coordinates("D:\\Experiment\\PLS\\Release\\checkins-10162.txt", coordinate_vector);
#else//相当于将用户中涉及的所有位置，都拿出来作为候选地的信息
		data_dump.load_string_checkin_without_time_to_coordinates("D:\\Experiment\\PLS\\Release\\checkins-2321.txt", coordinate_vector);
#endif	
		end_time = std::chrono::high_resolution_clock::now();
	//	std::cout << coordinate_vector.size() << std::endl;
		std::set<geo_coordinate> coordinate_set;
		std::vector<geo_coordinate>::const_iterator iter_venue = coordinate_vector.cbegin();
		for (; iter_venue != coordinate_vector.cend(); ++iter_venue)
		{
			coordinate_set.insert(*iter_venue);
		}

#if PICK_FROM_UNIQUE//利用了set数组关键字的唯一性，使得重复的位置没有重新写入coordinate_vector中
		coordinate_vector.clear();
		std::set<geo_coordinate>::const_iterator iter_unique_venue = coordinate_set.cbegin();
		for (; iter_unique_venue != coordinate_set.cend(); ++iter_unique_venue)
		{
			coordinate_vector.push_back(*iter_unique_venue);
		}
#endif

		std::cout << "2  checkins.txt has been loaded. " << coordinate_set.size() << " different venues." << std::endl;
		std::cout << "3  Loading took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << "ms" << std::endl;

		std::cout << "\n4  Generating " << candidate_count << " candidate locations." << std::endl;
		begin_time = std::chrono::high_resolution_clock::now();
		query_generator::random(coordinate_vector, candidate_vector, candidate_count);//从用户文件中读取的位置信息，然后在其中随机挑选candidate_count个位置
		end_time = std::chrono::high_resolution_clock::now();
		std::cout << "5  " << candidate_count << " candidate locations have been generated." << std::endl;
		std::cout << "6\tGenerating took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << "ms" << std::endl;

		std::cout << "\n7  Dumping candidate data to candidates.txt." << std::endl;
		begin_time = std::chrono::high_resolution_clock::now();

#if DO_50_TIMES
		data_dump.dump_string_coordinates(cand_file.c_str(), candidate_vector);
#else
		data_dump.dump_string_coordinates("D:\\Experiment\\PLS\\Release\\candidates.txt", candidate_vector);
#endif

		end_time = std::chrono::high_resolution_clock::now();
		std::cout << "8  candidates.txt has been dumped." << std::endl;
		std::cout << "9  Dumping took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << "ms" << std::endl;

		std::cout << "\n<<<<< CANDIDATES_GENERATION <<<<<\n" << std::endl;
#endif // CANDIDATES_GENERATION


#if CHECKINS_EXCLUDING
		std::cout << ">>>>> CHECKINS_EXCLUDING >>>>>\n" << std::endl;

		std::cout << "1  Loading and dumping check-in data excluding candidates." << std::endl;
		begin_time = std::chrono::high_resolution_clock::now();

#if GEN_FROM_10162
#if DO_50_TIMES
		std::string checkins_file = "D:\\Experiment\\PLS\\Release\\checkins-" + t_str + ".txt";
		data_dump.load_dump_checkins_excluding_candidates("D:\\Experiment\\PLS\\Release\\checkins-10162.txt", cand_file.c_str(), checkins_file.c_str());
#else
		data_dump.load_dump_checkins_excluding_candidates("D:\\Experiment\\PLS\\Release\\checkins-10162.txt", "D:\\Experiment\\PLS\\Release\\candidates.txt",
			"D:\\Experiment\\PLS\\Release\\checkins.txt");
#endif
#else
#if DO_50_TIMES
		std::string checkins_file = "D:\\Experiment\\PLS\\Release\\checkins-" + t_str + ".txt";
		data_dump.load_dump_checkins_excluding_candidates("D:\\Experiment\\PLS\\Release\\checkins-2321.txt", cand_file.c_str(), checkins_file.c_str());
#else
		data_dump.load_dump_checkins_excluding_candidates("D:\\Experiment\\PLS\\Release\\checkins-2321.txt", "D:\\Experiment\\PLS\\Release\\candidates.txt",
			"D:\\Experiment\\PLS\\Release\\checkins.txt");
#endif
#endif

		end_time = std::chrono::high_resolution_clock::now();
		std::cout << "2  Check-in data has been excluded by candidates." << std::endl;
		std::cout << "3  Loading and dumping took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << "ms" << std::endl;

		std::cout << "\n<<<<< CHECKINS_EXCLUDING <<<<<\n" << std::endl;
#endif // CHECKINS_EXCLUDING


#if DATA_LOADING
		std::cout << ">>>>> DATA_LOADING >>>>>\n" << std::endl;

		std::map<std::string, std::vector<geo_point>> checkin_map;
		std::map<std::string, std::vector<geo_point>> checkin_map_origin;

		std::cout << "1  Loading users' check-in data from checkins.txt to a map." << std::endl;
		begin_time = std::chrono::high_resolution_clock::now();

#if DO_50_TIMES
		data_dump.load_string_checkin_without_time_by_user(checkins_file.c_str(), checkin_map);
#else
		data_dump.load_string_checkin_without_time_by_user("D:\\Experiment\\PLS\\Release\\checkins.txt", checkin_map);
#endif

#if GEN_FROM_10162
		data_dump.load_string_checkin_without_time_by_user("D:\\Experiment\\PLS\\Release\\checkins-10162.txt", checkin_map_origin);
#else
		data_dump.load_string_checkin_without_time_by_user("D:\\Experiment\\PLS\\Release\\checkins-2321.txt", checkin_map_origin);
#endif
		end_time = std::chrono::high_resolution_clock::now();
		std::cout << "2  checkins.txt has been loaded." << std::endl;
		std::cout << "3  Loading took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << "ms" << std::endl;

		//std::cout << checkin_map.size() << std::endl;

		candidate_vector.clear();

		std::cout << "\n4  Loading candidate data from candidates.txt to a vector." << std::endl;
		begin_time = std::chrono::high_resolution_clock::now();

#if DO_50_TIMES
		data_dump.load_string_coordinates(cand_file.c_str(), candidate_vector);
#else
		data_dump.load_string_coordinates("D:\\Experiment\\PLS\\Release\\candidates.txt", candidate_vector);
#endif

		end_time = std::chrono::high_resolution_clock::now();
		std::cout << "5  candidates.txt has been loaded." << std::endl;
		std::cout << "6  Loading took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << "ms" << std::endl;

		std::cout << "\n<<<<< DATA_LOADING <<<<<\n" << std::endl;
#endif

		long double tau = 0.7;


#if DO_50_TIMES//
	}	// For times
#endif

				pino algorithm;
				algorithm.is_box_inf = true;
				int k = 10;//the number of target candidates.
				std::map<geo_coordinate, int> result;
				std::cout << "1  Setting data." << std::endl;
				//algorithm.set_probability_function(probability_functions::linear);
				//algorithm.set_probability_inverse_function(probability_functions::linear_inverse);
				algorithm.set_probability_function(probability_functions::power_law);
				algorithm.set_probability_inverse_function(probability_functions::power_law_inverse);
				algorithm.set_threshold(tau);
				algorithm.set_checkin_data(checkin_map); 			
				algorithm.set_candidate_data(candidate_vector);

				std::cout << "\n2  Executing NA algorithm." << std::endl;
				algorithm.algo_opt = pino::AT_NA;
				std::vector<geo_coordinate> result_na;
				unsigned max_inf_na = 0;
				begin_time = std::chrono::high_resolution_clock::now();
				algorithm.prepare();


				std::cout << "\n14  New PIN-VO algorithm." << std::endl;
				algorithm.algo_opt = pino::AT_NEW_VO;
				std::vector<geo_coordinate> result_new_vo;
				unsigned max_inf_new_vo = 0, IA_new_vo = 0, NIB_new_vo = 0;
				begin_time = std::chrono::high_resolution_clock::now();
				long long time_new_vo;
				algorithm.prepare();

				//Baseline Algorithm.
				std::vector<geo_coordinate> result_base;
				result_base.clear();
				std::cout << "\n22 Executing Baseline algorithm." << std::endl;
				begin_time = std::chrono::high_resolution_clock::now();
				std::cout << "k=" << k << "-----------" << "candidate=" << candidate_count << std::endl;
				algorithm.prepare();
				algorithm.Baseline(result_base, k);
				end_time = std::chrono::high_resolution_clock::now();
				time_new_vo = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();
				std::cout << "20 Executing Baseline algorithm took " << time_new_vo << "ms" << std::endl;
				int count_base;
				count_base = algorithm.test_traj1(result_base, k);
				std::cout << count_base << std::endl << std::endl;


				result.clear();
				std::cout << "\n17 Executing Top-k New PIN-VO algorithm." << std::endl; //GreedyP algorithm
				begin_time = std::chrono::high_resolution_clock::now();
				algorithm.prepare();
				int count_na_k = algorithm.pin_vo_k(result, k);
				end_time = std::chrono::high_resolution_clock::now();
				time_new_vo = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();
				std::cout << "18 Executing Top-k New PIN-VO took " << time_new_vo << "ms" << std::endl;


				
				for (int count = 0; count < 10; count++)
				{
					std::cout << "the count of result:" << count << std :: endl;
					result.clear();
					std::cout << "\n21 Executing num bit vector New PIN-VO algorithm." << std::endl;
					begin_time = std::chrono::high_resolution_clock::now();
					std::cout << "k=" << k << "------------" << "number of bit vector:" << HASH_NUM << "-----------" << "candidate=" << candidate_count << std::endl;
					algorithm.prepare();
					algorithm.pin_vo_num_vitVector(result, k);
					end_time = std::chrono::high_resolution_clock::now();
					time_new_vo = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();
					std::cout << "20 Executing num bit vector New PIN-VO took " << time_new_vo << "ms" << std::endl;
					int count_bit;
					count_bit = algorithm.test_traj(result, k); // calculate the number of objects influenced by target candidates.
					std::cout << count_bit << std::endl<<std::endl;

				}

	std::cout << "\nInput \'e\' and Enter to exit." << std::endl;
	char input = '\0';
	while (input != 'e')
		std::cin >> input;

	return 0;
}
//比较各种算法的运行时间，影响用户的最大值，IA，IB等信息，并且打印出每种算法选择的候选地
