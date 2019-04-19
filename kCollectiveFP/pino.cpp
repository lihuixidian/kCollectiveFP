#include "stdafx.h"
#include "pino.h"
#include<algorithm>
#include <queue>
#include <chrono>//��ʱ���йص�ͷ�ļ�
#include<cstdlib>
#include<ctime>
#define random(a,b) (rand()%(b-a+1)+a)



void pino::prepare()//�ж�Ҫ���е��㷨��Ȼ��ֱ�Ը��Ե����ݽṹ���г�ʼ������
{
	if (algo_opt == AT_NA)
		return;

	if (algo_opt == AT_CAND)//�����ѧ��˵��reviewer��������Բ���
	{
		index_users_rtree();
		return;
	}

	mmr_map.clear();
	// pruned_by_IA = pruned_by_NIB = 0 will be inited in each algo.

	if (is_box_inf)
		index_box_users();
	else
		index_ring_users();

	index_candidates();
}


// Calculate samples count points, corresponding minMaxRadius values, then prepare info (init inf & potential number, user condition, MBRs) for each candidate.
void pino::prepare_candidates(SP_TYPE type)
{
	if (type != SPT_LINEAR && type != SPT_COUNTS)
		type = SPT_LINEAR; // Set to default value.

	candidate_info_map.clear();
	count_mmr_map.clear();

	std::map<std::string, std::vector<geo_point>>::const_iterator iter_user; // User iterator.
	if (type == SPT_LINEAR)
	{
		//out_file << "SPT_LINEAR\n";
		iter_user = checkin_map.cbegin();
		unsigned min_count = 100000, max_count = 0;
		for (; iter_user != checkin_map.cend(); ++iter_user) // Get min & max counts of users.
		{
			unsigned size = iter_user->second.size(); // Positions count.�ҵ��û�λ�����������ֵ����Сֵ
			if (size > max_count)
				max_count = size;
			if (size < min_count)
				min_count = size;
		}

		unsigned step_count = (max_count - min_count) / (points_count * 2); // Half of positions count each section.
		for (unsigned i = 0; i < points_count; ++i) // Cal all sample points.
		{
			// Calculate mmr.
			unsigned n = min_count + (i * 2 + 1) * step_count;
			long double pro_n = 1.0 - pow(1.0 - threshold, 1.0 / n);
			long double mmr = probability_inverse_function(pro_n);
			count_mmr_map[n] = mmr;
			//out_file << n << ", " << mmr << "\n";
		}
	}
	else if (type == SPT_COUNTS)
	{
		//out_file << "SPT_COUNTS\n";
		std::map<unsigned, unsigned> user_count_map; // <positions count, user count>.
		iter_user = checkin_map.cbegin();
		for (; iter_user != checkin_map.cend(); ++iter_user) // Get the number of users for each count.
		{
			unsigned size = iter_user->second.size(); // Positions count.
			if (user_count_map.find(size) != user_count_map.end())//ͳ�Ʋ�ͬλ�ø������û�����
				user_count_map[size]++;
			else
				user_count_map[size] = 1;
		}

		unsigned user_count = checkin_map.size() / (points_count * 2); // Half of users count for each section.
		unsigned step_count = user_count;
		unsigned count_accumulation = 0, count_left = 0, count_right = 0;
		std::map<unsigned, unsigned>::const_iterator iter_count = user_count_map.cbegin();
		for (; iter_count != user_count_map.cend(); ++iter_count)
		{
			count_accumulation += iter_count->second;

			if (count_accumulation < step_count) // Not reach half of section.
				count_left = step_count - count_accumulation; // Be short of.
			else // Exceed (or exactly at) half of section.
			{
				count_right = count_accumulation - step_count; // Beyond (or zero).
				if (count_left < count_right) // Tend to less side (left).
				{
					std::map<unsigned, unsigned>::const_iterator iter_left = iter_count;
					--iter_left; // Previous count.

					// Calculate mmr.
					unsigned n = iter_left->first;
					long double pro_n = 1.0 - pow(1.0 - threshold, 1.0 / n);
					long double mmr = probability_inverse_function(pro_n);
					count_mmr_map[n] = mmr;
					//out_file << n << ", " << mmr << "\n";
				}
				else // Tend to less side (right or zero).
				{
					// Calculate mmr.
					unsigned n = iter_count->first; // Current count.
					long double pro_n = 1.0 - pow(1.0 - threshold, 1.0 / n);
					long double mmr = probability_inverse_function(pro_n);
					count_mmr_map[n] = mmr;
					//out_file << n << ", " << mmr << "\n";
				}				

				step_count += user_count * 2; // Update to next section.
				count_left = step_count - count_accumulation; // New be short of.
				count_right = 0;
			}
		}
	}

	// Init candidate info for condition checking and max heap.
	unsigned user_count = checkin_map.size(); // Init potential inf value.maxInf
	std::vector<geo_coordinate>::const_iterator iter_cand = candidate_vector.cbegin();
	for (; iter_cand != candidate_vector.cend(); ++iter_cand) // Iterate each candidate.����ÿһ����ѡ�ص�
	{
		candidate_info_ptr cand_ptr(new candidate_info()); // Construct a candidate_info object shared_pointer.

		// Init candidate info for inf & potential values.
		cand_ptr->inf = 0;
		cand_ptr->potential = user_count;

		// Init candidate info for condition checking.
		for (iter_user = checkin_map.cbegin(); iter_user != checkin_map.cend(); ++iter_user) // Iterate each user.
			cand_ptr->condition_map[iter_user->first] = 0; // Init each user condition.

		// Init candidate info for circle MBRs.�Ժ�ѡ��ΪԲ�ģ�minMaxRadiusΪ�뾶��Բ���ڽӾ��Σ���������������
		std::map<unsigned, long double>::const_iterator iter_mmr = count_mmr_map.cbegin();
		for (; iter_mmr != count_mmr_map.cend(); ++iter_mmr) // Iterate each sample points (and the corresponding mmr).һ���Ժ�ѡ��ΪԲ�ģ���ͬ��minMaxRadius��n��t��Ϊ�뾶���ڽӾ���
		{
			geo_box mbr;
			mbr.max_corner() = geo_offset(geo_offset(iter_cand->get_coordinate(), iter_mmr->second, northward), iter_mmr->second, eastward);
			mbr.min_corner() = geo_offset(geo_offset(iter_cand->get_coordinate(), iter_mmr->second, southward), iter_mmr->second, westward);
			cand_ptr->circle_mbr_map[iter_mmr->first] = mbr; // Init each circle MBR.ÿ���û���λ����N��ͬ�����Ժ�ѡ�ص��ڽӾ����кö��
		}

		candidate_info_map[*iter_cand] = cand_ptr; // Associate candidate and its info shared_pointer.
	}
}


void pino::execute_algorithm(std::vector<geo_coordinate>& query_result, unsigned& max_inf)//�����㷨
{
	switch (algo_opt)
	{
	case AT_NA://����NA�㷨navie
		na(query_result, max_inf);
		break;
	case AT_VO://����PIN_VO*�㷨û�м�֦����
		vo(query_result, max_inf);
		break;
	case AT_PIN://����PIN�㷨2
		if (is_box_inf)
			pin_box(query_result, max_inf);        
		else
			;
		break;
	case AT_PIN_VO://�ɵ�PIN-VO�㷨3
		if (is_box_inf)
			pin_vo_box(query_result, max_inf);
		else
			;
		break;
	case AT_NEW_VO://�µ�PIN-VO�㷨3
		if (is_box_inf)
			new_vo_box(query_result, max_inf);
		else
			;
		break;
	case AT_CAND://�������
		cand(query_result, max_inf);
		break;
	default:
		break;
	}
}


// <PS> Assume all input data have been set, therefore no data checking is implemented.
//���������û����ѡ�ص�Prc(Oi)��ѡ��Ӱ���û������Ǹ���ѡ����Ϊ���λ��
void pino::na(std::vector<geo_coordinate>& query_result, unsigned& max_inf)
{
	query_result.clear();
	max_inf = 0;

	// Init inf value for each candidate.
	std::map<geo_coordinate, unsigned> candidate_inf_map;
	std::vector<geo_coordinate>::const_iterator iter_candidate = candidate_vector.cbegin();
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate)
	{
		candidate_inf_map[*iter_candidate] = 0;//��ÿ����ѡ��Ӱ����û�������ʼ��Ϊ0
	}

	iter_candidate = candidate_vector.cbegin();
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate) // Iterate each candidate to calculate its inf value.
	{
		std::map<std::string, std::vector<geo_point>>::const_iterator iter_checkin = checkin_map.cbegin();
		for (; iter_checkin != checkin_map.cend(); ++iter_checkin) // Iterate each user.
		{
			long double product = 1.0; // "1 - product >= threshold".
			std::vector<geo_point>::const_iterator iter_coordinate = iter_checkin->second.cbegin();
			for (; iter_coordinate != iter_checkin->second.cend(); ++iter_coordinate) // Iterate each position of a user.
			{
				long double distance = geo_distance(iter_candidate->get_coordinate(), *iter_coordinate);
				product *= 1.0 - probability_function(distance);//����Prc��p1����Prc��p2����...��Prc��pn��
			}
			if (1.0 - product >= threshold)//Prc��O��
				++(candidate_inf_map[*iter_candidate]);
		} // Now, all users are checked for this candidate.

		// Update global max inf & current optimal candidate(s).
		if (candidate_inf_map[*iter_candidate] > max_inf)//�����ֱȵ�ǰmax_inf��������ĺ�ѡ�أ����query_result����������ѹ���µĺ�ѡ��
		{
			max_inf = candidate_inf_map[*iter_candidate];
			query_result.clear();
			query_result.push_back(*iter_candidate);
		}
		else if (candidate_inf_map[*iter_candidate] == max_inf)//����ǰ��ѡ����Ӱ����û������뵱ǰ���ֵ��ͬ���򽫷�����ѡ��ѹ��
			query_result.push_back(*iter_candidate);//����ǰ��ѡ��Ӱ���û��������ڵ�ǰ���ֵ����ʲô�������������鿴�µĺ�ѡ��
	}
	int a = 1;
}


// Remark: C_POTENTIAL & C_TO_VERIFY are not used in this function.�㷨����ʵ��
//������㷨�У�û�ж�cadidate_tuple_map��candidate_tuple�еĵ�����Ԫ��Ҳ������Ҫ�ú�ѡ����֤��vector�û����г�ʼ��
//�ڱ���ÿһ���û�ʱ���ҵ��ڸ��û�IA��NIB�ڵĺ�ѡ���Ժ�ֱ��ֻ��������minInf��maxInf��ֵ��
void pino::pin_box(std::vector<geo_coordinate>& query_result, unsigned& max_inf)//query_result������������Ѻ�ѡλ�õ�������max_inf�Ǻ�ѡ����Ӱ��������û�����
{
	query_result.clear();
	max_inf = 0;

	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size();//��ѡ�ص�����

	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.����ÿ���û�
	{
		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // Size of candidate in IA.
		bool inf_valid = iter_user->second.get<U_INF_VALID>();
		if (inf_valid)
		{
			// Query candidates inside IA.
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),//query���ص��Ǻ�ѡ�����û���IA��Χ�ڵĸ���
				std::inserter(inf_candidates, inf_candidates.begin()));//�ڸô�inf_candidates�д�ŵľ��Ǻ�ѡ����IA�еĵص�
			pruned_by_IA += inf_size; // Accumulate inf values.

			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				unsigned inf = ++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment inf value by 1 for each candidate inside IA.

				// Update max candidate.
				if (inf > max_inf)
				{
					max_inf = inf;
					query_result.clear();
					query_result.push_back(*iter_inf);
				}
				else if (inf == max_inf)
					query_result.push_back(*iter_inf);
			}
		}

		// Query potential candidates inside NIB.
		std::vector<geo_coordinate> potential_candidates;//������¼�û�NIB��Χ���ڵĺ�ѡ�ؼ���
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		//candidate_rtree.query()�������ã�����NIB�߽緶Χ�ڵĺ�ѡ�ظ������أ����ҽ���ѡ����Ϣ������Potential_candidates��
		pruned_by_NIB += candidates_count - potential_size;; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{
			if (inf_size > 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
				continue;
			//������������IA��NIB���ڵĺ�ѡ��
			candidate_tuple* candidate_potential_tuple = &(candidate_tuple_map[potential_candidates[i]]); // The tuple of potential candidate to be checked.

			geo_point candidate(potential_candidates[i].get_coordinate()); // Potential candidate coordinate.

			long double product = 1.0; // "1 - product >= threshold".
			std::vector<geo_point>::const_iterator iter_coordinate = iter_user->second.get<U_POINTS_PTR>()->cbegin();
			for (; iter_coordinate != iter_user->second.get<U_POINTS_PTR>()->cend(); ++iter_coordinate) // Iterate each position of the user.
			{
				long double distance = geo_distance(candidate, *iter_coordinate);
				product *= 1.0 - probability_function(distance);
			}

			if (1.0 - product >= threshold)
			{
				unsigned cur_inf = ++(candidate_potential_tuple->get<C_INF>());
				if (cur_inf > max_inf)
				{
					max_inf = cur_inf;
					query_result.clear();
					query_result.push_back(potential_candidates[i]);
				}
				else if (cur_inf == max_inf)
					query_result.push_back(potential_candidates[i]);
			}
		}
	}
}


// Remark: C_TO_VERIFY is not used in this function.�ú���ʵ���㷨������˵�Ǿɰ�
void pino::pin_vo_box(std::vector<geo_coordinate>& query_result, unsigned& max_inf)
{
	query_result.clear();
	max_inf = 0; // maxminInf.

	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size();

	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		// Decrease "maxInf" by 1 for all candidates.
		std::map<geo_coordinate, candidate_tuple>::iterator iter = candidate_tuple_map.begin();
		for (; iter != candidate_tuple_map.end(); ++iter)
			--(iter->second.get<C_POTENTIAL>());

		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // Size of candidate in IA.
		bool inf_valid = iter_user->second.get<U_INF_VALID>();
		if (inf_valid)
		{
			// Query candidates inside IA.
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),
				std::inserter(inf_candidates, inf_candidates.begin()));
			pruned_by_IA += inf_size; // Accumulate inf values.

			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				unsigned inf = ++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment "minInf" by 1 for each candidate inside IA.

				// Update max candidate.
				if (inf > max_inf)
				{
					max_inf = inf;
					query_result.clear();
					query_result.push_back(*iter_inf);
				}
				else if (inf == max_inf)
					query_result.push_back(*iter_inf);
			}
		}

		// Query potential candidates inside NIB.�൱����NIB֮��ĺ�ѡ��û�����κδ���
		std::vector<geo_coordinate> potential_candidates;
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		pruned_by_NIB += candidates_count - potential_size;; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{		
			// Upper and lower bounds.
			candidate_tuple* candidate_potential_tuple = &(candidate_tuple_map[potential_candidates[i]]); // The tuple of potential candidate to be checked.
			if (candidate_potential_tuple->get<C_POTENTIAL>() + 1 <= max_inf) // Candidate potential value <= global max inf, then no need to increment it by 1.
				continue;

			if (inf_size > 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
			{
				++(candidate_potential_tuple->get<C_POTENTIAL>()); // This candidate is in IA, then we bring back 1 to potential value.ǰ�涼������1�������ֶ���1
				continue;
			}

			geo_point candidate(potential_candidates[i].get_coordinate()); // Potential candidate coordinate, which is between IA and NIB.

			long double product = 1.0; // "1 - product >= threshold".
			std::vector<geo_point>::const_iterator iter_coordinate = iter_user->second.get<U_POINTS_PTR>()->cbegin();
			for (; iter_coordinate != iter_user->second.get<U_POINTS_PTR>()->cend(); ++iter_coordinate) // Iterate each position of the use.
			{
				long double distance = geo_distance(candidate, *iter_coordinate);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold)	// Early stopping.
				{
					unsigned cur_inf = ++(candidate_potential_tuple->get<C_INF>()); // Whenever "1.0 - product >= threshold" holds, the candidate must inf the user.
					++(candidate_potential_tuple->get<C_POTENTIAL>()); // This candidate inf the user, then we bring back 1 to potential value.
					
					if (cur_inf > max_inf) // Update global (max_inf) maxminInf.
					{
						max_inf = cur_inf;
						query_result.clear();
						query_result.push_back(potential_candidates[i]);
					}
					else if (cur_inf == max_inf)
						query_result.push_back(potential_candidates[i]);

					break; // No need to check remainder positions of the user.�������2����ǰ����ѭ��
				}
			}				
		} // End of potential candidates iteration.
	}
}

//�㷨3���°棬Ҳ���������е��㷨��
void pino::new_vo_box(std::vector<geo_coordinate>& query_result, unsigned& max_inf)
{
	query_result.clear();
	max_inf = 0; // maxminInf.

	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size(); // Candidates count.

	//auto begin_time = std::chrono::high_resolution_clock::now();
	//auto end_time = std::chrono::high_resolution_clock::now();

	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // The number of candidates inside IA.
		if (iter_user->second.get<U_INF_VALID>()) // Check IA is valid or not.
		{
			// Query candidates inside IA.
			//begin_time = std::chrono::high_resolution_clock::now();
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),
				std::inserter(inf_candidates, inf_candidates.begin()));
			//end_time = std::chrono::high_resolution_clock::now();
			pruned_by_IA += inf_size; // Accumulate inf values.

			// Iterate each candidate in IA.
			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment inf value by 1 for each candidate inside IA.
			}//������IA��Χ�ڵĺ�ѡ�أ�ֻ��Ҫ��minInf+1����
		}

		// Query potential candidates inside NIB.
		std::vector<geo_coordinate> potential_candidates;
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		pruned_by_NIB += candidates_count - potential_size; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{
			if (inf_size != 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
				continue;

			candidate_tuple_map[potential_candidates[i]].get<C_TO_VERIFY>().push_back(iter_user->first); // Append this user as a potential user for the candidate.
		}//candidate_tuple_map�����ѡ����Ҫ��֤�û�����Ϣ
	}
	//long long time_cand_counts = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();
	//out_file << "rtree time=" << time_cand_counts << "ns\n";

	// Use a max_heap to rank candidates by potential inf and min inf.
	std::priority_queue<CandidateIter> max_heap;
	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		unsigned p_value = iter_tuple->second.get<C_POTENTIAL>() = iter_tuple->second.get<C_INF>() + iter_tuple->second.get<C_TO_VERIFY>().size();	// Potential inf.
		max_heap.push(CandidateIter(p_value, iter_tuple->second.get<C_INF>(), iter_tuple));//�Զѽ��г�ʼ�������Ҹ�����maxInf��ֵ
	}

	//out_file << "New PIN-VO\n";
	//begin_time = std::chrono::high_resolution_clock::now();
	// Iterate candidates based on (potential and min infs) max heap.����Ҫ��֤�ĺ�ѡ�ؽ��д���
	while (!max_heap.empty())
	{
		// Upper and lower bounds.
		if (max_heap.top().potential_value < max_inf)
			break;
		//out_file << max_heap.top().iter_map->first.longitude() << ", " << max_heap.top().iter_map->first.latitude() << ", g=" << max_inf << ", p=" << max_heap.top().potential_value << ", i=" << max_heap.top().min_inf_value << "\n";

		geo_point candidate = max_heap.top().iter_map->first.get_coordinate(); // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = max_heap.top().iter_map->second.get<C_TO_VERIFY>().cbegin(),
												 iter_end = max_heap.top().iter_map->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
												   iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate, *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(max_heap.top().iter_map->second.get<C_INF>()); // Update "minInf" value.
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}
			
			if (!is_inf // Not inf this user.
				&& --(max_heap.top().iter_map->second.get<C_POTENTIAL>()) < max_inf) // Potential <= global max inf.
				break; // The candidate is not optimal.�����ѡ�ز���Ӱ����û�����maxInf-1�����������1����ֹͣ�Ըú�ѡ�ص���֤
		}

		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);

		max_heap.pop(); // Pop next candidate to heap top.
	}
	//end_time = std::chrono::high_resolution_clock::now();
	//time_cand_counts = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();
	//out_file << "time=" << time_cand_counts << "\n";
}


// Remark: C_TO_VERIFY are not used in this function.ֻ�������Ż����ԣ�û�У����ü�֦����
void pino::vo(std::vector<geo_coordinate>& query_result, unsigned& max_inf)
{
	query_result.clear();
	max_inf = 0; // maxminInf.

	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		std::map<geo_coordinate, candidate_tuple>::iterator iter_cand = candidate_tuple_map.begin();
		for (; iter_cand != candidate_tuple_map.end(); ++iter_cand) // Iterate each candidate.
		{
			// Upper and lower bounds.
			if (iter_cand->second.get<C_POTENTIAL>() <= max_inf) // Candidate potential value <= global max inf, then no need to check it.
				continue;//ֱ�ӱ�����һ����ѡ��

			geo_point candidate(iter_cand->first.get_coordinate()); // Candidate coordinate.

			long double product = 1.0; // "1 - product >= threshold".
			bool is_inf = false;
			std::vector<geo_point>::const_iterator iter_coordinate = iter_user->second.get<U_POINTS_PTR>()->cbegin();
			for (; iter_coordinate != iter_user->second.get<U_POINTS_PTR>()->cend(); ++iter_coordinate) // Iterate each position of the use.
			{
				long double distance = geo_distance(candidate, *iter_coordinate);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold)	// Early stopping.
				{
					is_inf = true;
					unsigned cur_inf = ++(iter_cand->second.get<C_INF>()); // Whenever "1.0 - product >= threshold" holds, the candidate must inf the user.

					if (cur_inf > max_inf) // Update global (max_inf) maxminInf.
					{
						max_inf = cur_inf;
						query_result.clear();
						query_result.push_back(iter_cand->first);
					}
					else if (cur_inf == max_inf)
						query_result.push_back(iter_cand->first);

					break; // No need to check remainder positions of the user.
				}
			}

			if (!is_inf
				&& --(iter_cand->second.get<C_POTENTIAL>()) <= max_inf) // Decrease "maxInf" by 1, as it cannot inf the user.
				break;

		} // End of candidates iteration.
	}
}


void pino::index_box_users()
{
	user_box_tuple_map.clear();

	std::map<std::string, std::vector<geo_point>>::const_iterator iter_user = checkin_map.cbegin();//���ౣ���û�λ����Ϣ��map
	for (; iter_user != checkin_map.cend(); ++iter_user) // Iterate each user.
	{
		// xy multi points (positions of a user).
		geo_multi_point_ptr points_ptr(new geo_multi_point(iter_user->second.begin(), iter_user->second.end()));//ÿ���û�������λ����һ��multipoint��ָ��ָ��

		// Cardinality of positions.
		unsigned n = iter_user->second.size();//ÿ���û����е�λ������

		// minMaxRadius(T,n).
		double mmr = 0.0;
		std::map<unsigned, long double>::iterator iter_mmrs = mmr_map.find(n);//�ӹ�ϣ������λ�ø���Ϊn��minMaxRadius��ֵ
		if (iter_mmrs != mmr_map.end()) // Find a mmr that has been computed before.
		{
			mmr = iter_mmrs->second;
		}
		else // Cannot find, so compute it.�ڹ�ϣ����û���ҵ���Ӧn��minMaxRadius��ֵ�������¼���
		{
			long double pro_n = 1.0 - pow(1.0 - threshold, 1.0 / n);
			mmr = probability_inverse_function(pro_n);
			mmr_map[n] = mmr;//map���������鸳ֵ��ʽ
		}

		// MBR of the user.
		geo_box mbr = boost::geometry::return_envelope<geo_box>(*(points_ptr.get()));//��ĳһ���û�������λ�����ҵ�MBR

		// Inf bounds are valid or not.
		geo_point centroid;
		boost::geometry::centroid(mbr, centroid);//boost::geometry::centroid�������㼸�����ģ���һ������Ϊ������״�����û���MBR���ڶ���������Ϊ���ĵ�ֵ
		bool is_inf_bounds_valid = geo_distance(mbr.max_corner(), centroid) <= mmr;//���ĵ�mbr��ǵľ���С��mmr������������������������

		// The bound values are approximation viewing as a plane distance.
		geo_box inf_bound; // Inf bound.IA����
		if (is_inf_bounds_valid)
		{
			long double half_x = geo_distance(mbr.min_corner(), geo_point(centroid.get<0>(), mbr.min_corner().get<1>())); // Half x of MBR.
			long double half_y = geo_distance(mbr.min_corner(), geo_point(mbr.min_corner().get<0>(), centroid.get<1>())); // Half y of MBR.
			long double gamma = (acos(half_x / mmr) - asin(half_y / mmr)) / 2 + asin(half_y / mmr); // Angle at 1/8 arc point.
			long double gamma_x = cos(gamma) * mmr;//�ڽӾ��εĳ�
			long double gamma_y = sin(gamma) * mmr;//�ڽӾ��εĿ�
			inf_bound.max_corner() = geo_offset(geo_offset(mbr.min_corner(), gamma_y, northward), gamma_x, eastward); // Is anti-corner of min_corner.
			inf_bound.min_corner() = geo_offset(geo_offset(mbr.max_corner(), gamma_y, southward), gamma_x, westward); // Is anti-corner of max_corner.
		}
		geo_box potential; // Potential bound.NIB����
		potential.min_corner() = geo_offset(geo_offset(mbr.min_corner(), mmr, southward), mmr, westward);
		potential.max_corner() = geo_offset(geo_offset(mbr.max_corner(), mmr, northward), mmr, eastward);

		user_box_tuple_map[iter_user->first] = boost::make_tuple(n, points_ptr, is_inf_bounds_valid, inf_bound, potential);
	}
}


void pino::index_ring_users()
{
	user_ring_tuple_map.clear();

	std::map<std::string, std::vector<geo_point>>::const_iterator iter_user = checkin_map.cbegin();
	for (; iter_user != checkin_map.cend(); ++iter_user) // Iterate each user.
	{
		// xy multi points (positions of a user).
		geo_multi_point_ptr points_ptr(new geo_multi_point(iter_user->second.begin(), iter_user->second.end()));

		// Cardinality of positions.
		unsigned n = iter_user->second.size();

		// minMaxRadius(T,n).
		double mmr = 0.0;
		std::map<unsigned, long double>::iterator iter_mmrs = mmr_map.find(n);
		if (iter_mmrs != mmr_map.end()) // Find a mmr that has been computed before.
		{
			mmr = iter_mmrs->second;
		}
		else // Cannot find, so compute it.
		{
			long double pro_n = 1.0 - pow(1.0 - threshold, 1.0 / n);
			mmr = probability_inverse_function(pro_n);
			mmr_map[n] = mmr;
		}

		// MBR of the user.
		geo_box mbr = boost::geometry::return_envelope<geo_box>(*(points_ptr.get()));

		// Inf bounds are valid or not.
		geo_point centroid;
		boost::geometry::centroid(mbr, centroid);
		bool is_inf_bounds_valid = geo_distance(mbr.max_corner(), centroid) <= mmr;

		// The bound values are approximation viewing as a plane distance.
		geo_ring inf_bound; // Inf bound.
		if (is_inf_bounds_valid)//���������ʽ�����ڽӶ����������IA
		{
			long double half_x = geo_distance(mbr.min_corner(), geo_point(centroid.get<0>(), mbr.min_corner().get<1>())); // Half x of MBR.
			long double half_y = geo_distance(mbr.min_corner(), geo_point(mbr.min_corner().get<0>(), centroid.get<1>())); // Half y of MBR.
			long double alpha = acos(half_x / mmr);
			long double alpha_y = mmr * sin(alpha);
			long double beta = asin(half_y / mmr);
			long double beta_x = mmr * cos(beta);
			long double gamma = (alpha - beta) / 2 + beta; // Angle at 1/8 arc point.
			long double gamma_x = cos(gamma) * mmr;
			long double gamma_y = sin(gamma) * mmr;

			// Compute x/y for 8 points.
			geo_point clock12 = geo_offset(geo_point(centroid.get<0>(), mbr.min_corner().get<1>()), alpha_y, northward);
			geo_point clock1 = geo_offset(geo_offset(mbr.min_corner(), gamma_y, northward), gamma_x, eastward); // Is anti-corner of min_corner.
			geo_point clock3 = geo_offset(geo_point(mbr.min_corner().get<0>(), centroid.get<1>()), beta_x, eastward);
			geo_point clock5 = geo_offset(geo_offset(geo_point(mbr.min_corner().get<0>(), mbr.max_corner().get<0>()), gamma_y, southward), gamma_x, eastward);
			geo_point clock6 = geo_offset(geo_point(centroid.get<0>(), mbr.max_corner().get<1>()), alpha_y, southward);
			geo_point clock7 = geo_offset(geo_offset(mbr.max_corner(), gamma_y, southward), gamma_x, westward); // Is anti-corner of max_corner.
			geo_point clock9 = geo_offset(geo_point(mbr.max_corner().get<0>(), centroid.get<1>()), beta_x, westward);
			geo_point clock11 = geo_offset(geo_offset(geo_point(mbr.max_corner().get<0>(), mbr.min_corner().get<1>()), gamma_y, northward), gamma_x, westward);
			
			// Clockwise and closed.
			boost::geometry::append(inf_bound, clock12);
			boost::geometry::append(inf_bound, clock1);
			boost::geometry::append(inf_bound, clock3);
			boost::geometry::append(inf_bound, clock5);
			boost::geometry::append(inf_bound, clock6);
			boost::geometry::append(inf_bound, clock7);
			boost::geometry::append(inf_bound, clock9);
			boost::geometry::append(inf_bound, clock11);
			boost::geometry::append(inf_bound, clock12);
		}
		geo_box potential; // Potential bound��Ȼ���Ǹ���Ӿ���
		potential.min_corner() = geo_offset(geo_offset(mbr.min_corner(), mmr, southward), mmr, westward);
		potential.max_corner() = geo_offset(geo_offset(mbr.max_corner(), mmr, northward), mmr, eastward);

		user_ring_tuple_map[iter_user->first] = boost::make_tuple(n, points_ptr, is_inf_bounds_valid, inf_bound, potential);
	}
}


void pino::index_candidates()//��ʼ����ѡ�ص�tuple��ͺ�ѡ��R��
{
	candidate_rtree.clear();
	candidate_tuple_map.clear();

	unsigned user_count = checkin_map.size(); // Init potential inf value.

	std::vector<geo_coordinate>::const_iterator iter_candidate = candidate_vector.cbegin();//����̳е��Ĺ��ں�ѡ�ص���Ϣ
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate)
	{
		candidate_rtree.insert(iter_candidate->get_coordinate()); // Construct candidate R*-tree.������R��
		candidate_tuple_map[*iter_candidate] = boost::make_tuple(0, user_count, std::vector<std::string>()); // Init candidate inf values and vector.
		//��ʾһ���յ�vector����std::vector<std::string>()
	}
}


void pino::cand(std::vector<geo_coordinate>& query_result, unsigned& max_inf)//����Theorem1��Theorem2�����жԺ�ѡ�ص�ɸѡ
{
	query_result.clear();
	max_inf = 0; // maxminInf.

	std::priority_queue<CandidateInfoIter> max_heap; // Use a max_heap to rank candidates by potential inf and min inf.

	//auto begin_time = std::chrono::high_resolution_clock::now();
	//auto end_time = std::chrono::high_resolution_clock::now();

	std::map<geo_coordinate, candidate_info_ptr>::const_iterator iter_cand_info = candidate_info_map.cbegin(),
																 iter_cand_info_end = candidate_info_map.cend();
	for (; iter_cand_info != iter_cand_info_end; ++iter_cand_info) // Iterate each candidate info.
	{
		std::map<unsigned, geo_box>::const_iterator iter_mbr = iter_cand_info->second->circle_mbr_map.cbegin(),
													iter_mbr_end = iter_cand_info->second->circle_mbr_map.cend();
		for (; iter_mbr != iter_mbr_end; ++iter_mbr) // Iterate each circle MBR by count point for each candidate info.������ѡ�ز�ͬ�뾶��map
		{
			unsigned cur_count = iter_mbr->first; // Current count value (n).

			// Query users intersects c.MBR(n).
			std::vector<user_rtree_tuple> intersects_users;
			//begin_time = std::chrono::high_resolution_clock::now();
			user_rtree.query(boost::geometry::index::intersects(iter_mbr->second), std::back_inserter(intersects_users));
			//�ҵ����û�MBR�ڵ�ǰ��ѡ�ص�MBR�е������û�
			//end_time = std::chrono::high_resolution_clock::now();

			std::set<std::string> intersects_users_id; // In order to find user id faster.//�������ں�ѡ��ΪԲ�ĵ�minMaxRadiusΪ�뾶��Բ�е��û���Ϣ���и���
			for (unsigned i = 0; i < intersects_users.size(); ++i) // Iterate each user intersects c.MBR(n).
			{
				std::string user_id = intersects_users[i].get<UR_USER>();//User id
				intersects_users_id.insert(user_id); // Init user id set.
				if (iter_cand_info->second->condition_map[user_id] == 0) // The user still needs to be verified (condition = 0).
				{
					if (intersects_users[i].get<UR_NUM>() >= cur_count) // Actual number m >= n, must inf.�������n����minMaxRadius����
					{
						// maxDist vs. mmr checking.
						bool le_than = max_dist_le_than_mmr(iter_cand_info->first.get_coordinate(), intersects_users[i].get<UR_MBR>(), cur_count);
						if (le_than) // maxDist <= mmr, means inf.
						{
							iter_cand_info->second->condition_map[user_id] = 1; // The user is influenced by iter_cand_info (condition = 1).
							++(iter_cand_info->second->inf); // Inf value increments by 1.
						}
					} // For "else" (m < n), is still potential (condition = 0).
				} // For "else" (condition = 1 or -1).
			}

			// Check users not intersects c.MBR(n).Theorem2
			std::map<std::string, std::vector<geo_point>>::const_iterator iter_user = checkin_map.cbegin(); // User iterator.
			for (; iter_user != checkin_map.cend(); ++iter_user) // Iterate each user.
			{
				std::string user_id = iter_user->first;
				if (intersects_users_id.find(user_id) == intersects_users_id.end()) // The user is not intersects c.MBR(n).
				{
					if (iter_cand_info->second->condition_map[user_id] == 0) // The user still needs to be verified (condition = 0).
					{
						if (iter_user->second.size() <= cur_count) // Actual number m <= n, must not inf.����������û̫��
						{
							iter_cand_info->second->condition_map[user_id] = -1; // The user is not influenced by iter_cand_info (condition = -1).
							--(iter_cand_info->second->potential); // Potential value decrements by 1.
						} // For "else" (m > n), is still potential (condition = 0).
					} // For "else" (condition = 1 or -1).
				}
			}

		}  // End of iteration of each circle MBR.

		// Push each candidate into max_heap.
		max_heap.push(CandidateInfoIter(iter_cand_info->second->potential, iter_cand_info->second->inf, iter_cand_info));
	}
	//long long time_cand_counts = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();
	//out_file << "rtree time=" << time_cand_counts << "ns\n";
	
	//begin_time = std::chrono::high_resolution_clock::now();
	// Iterate candidates based on (potential and min infs) max_heap.
	while (!max_heap.empty())
	{
		// Upper and lower bounds.
		if (max_heap.top().potential_value < max_inf)
			break;
		//out_file << "g=" << max_inf << ", p=" << max_heap.top().potential_value << ", i=" << max_heap.top().min_inf_value << "\n";

		geo_point candidate = max_heap.top().iter_map->first.get_coordinate(); // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::map<std::string, int>::const_iterator iter_user = max_heap.top().iter_map->second->condition_map.cbegin(),
												   iter_end = max_heap.top().iter_map->second->condition_map.cend();
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.//�ҵ���Ȼ��Ҫ��֤���û�
		{
			if (iter_user->second != 0) // This user is determinate before.
				continue;

			std::vector<geo_point>::const_iterator iter_pos = checkin_map[iter_user->first].cbegin(),
												   iter_pos_end = checkin_map[iter_user->first].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate, *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(max_heap.top().iter_map->second->inf); // Update "minInf" value.
					break; // Jump out and compute next user (next iter_user).
				}				
			}

			if (!is_inf // Not inf this user.
				&& --(max_heap.top().iter_map->second->potential) < max_inf) // Potential <= global max inf.
				break; // The candidate is not optimal.
			//���˵��ǰ��ѡ����Ӱ���û������ֵ<maxInf������ڴ������˵��ʣ�¶��е�Ԫ��û�б�Ҫ�ٽ�����֤
		}

		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);

		max_heap.pop(); // Pop next candidate to heap top.
	}
	//end_time = std::chrono::high_resolution_clock::now();
	//time_cand_counts = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();
	//out_file << "time=" << time_cand_counts << "\n";
	//out_file.flush();
}


bool pino::max_dist_le_than_mmr(const geo_point& candidate, const geo_box& user, unsigned count)//�жϺ�ѡ�ؾ����û���Զ��ľ����Ƿ�<=minMaxRadius
{
	geo_point user_centroid, max_dist_corner;
	boost::geometry::centroid(user, user_centroid);

	if (candidate.get<0>() >= user_centroid.get<0>()) // c.lon >= u.lon.
		max_dist_corner.set<0>(user.min_corner().get<0>());
	else
		max_dist_corner.set<0>(user.max_corner().get<0>());
	if (candidate.get<1>() >= user_centroid.get<1>()) // c.lat >= u.lat.
		max_dist_corner.set<1>(user.min_corner().get<1>());
	else
		max_dist_corner.set<1>(user.max_corner().get<1>());

	return geo_distance(candidate, max_dist_corner) <= count_mmr_map[count];
}


void pino::index_users_rtree()
{
	user_rtree.clear();

	std::map<std::string, std::vector<geo_point>>::const_iterator iter_user = checkin_map.cbegin();
	for (; iter_user != checkin_map.cend(); ++iter_user) // Get min & max counts of users.
	{
		geo_multi_point points(iter_user->second.begin(), iter_user->second.end()); // xy multi points (positions of a user).
		unsigned n = iter_user->second.size(); // Cardinality of positions.
		geo_box mbr = boost::geometry::return_envelope<geo_box>(points); // MBR of the user.

		user_rtree.insert(boost::make_tuple(mbr, iter_user->first, n));	// Construct candidate R*-tree.
	}
}

void pino::na_k(std::map<geo_coordinate, int> & result, int k)
{	
	
	unsigned max_inf = 0;
	std::vector<geo_coordinate> query_result;
	// Init inf value for each candidate.
	std::map<geo_coordinate, unsigned> candidate_inf_map;
	std::vector<std::string> vec_info;//������Ӱ��ĵ�
	std::vector<geo_coordinate>::const_iterator iter_candidate = candidate_vector.cbegin();
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate)
	{
		candidate_inf_map[*iter_candidate] = 0;//��ÿ����ѡ��Ӱ����û�������ʼ��Ϊ0
	}
	
	iter_candidate = candidate_vector.cbegin();//������ѡ��vector
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate) // Iterate each candidate to calculate its inf value.
	{
		std::map<std::string, std::vector<geo_point>>::const_iterator iter_checkin = checkin_map.cbegin();
		
		vec_info.clear();
		for (; iter_checkin != checkin_map.cend(); ++iter_checkin) // Iterate each user.
		{
			long double product = 1.0; // "1 - product >= threshold".
			std::vector<geo_point>::const_iterator iter_coordinate = iter_checkin->second.cbegin();
			for (; iter_coordinate != iter_checkin->second.cend(); ++iter_coordinate) // Iterate each position of a user.
			{
				long double distance = geo_distance(iter_candidate->get_coordinate(), *iter_coordinate);
				product *= 1.0 - probability_function(distance);//����Prc��p1����Prc��p2����...��Prc��pn��
			}
			if (1.0 - product >= threshold)//Prc��O��
			{
				
				++(candidate_inf_map[*iter_candidate]);
			//	std::cout << iter_checkin->first << std::endl;
				vec_info.push_back(iter_checkin->first);
			}
		} // Now, all users are checked for this candidate.
		if (vec_info.size())
		{
			cand_info[*iter_candidate] = vec_info;

		}
		// Update global max inf & current optimal candidate(s).
		if (candidate_inf_map[*iter_candidate] > max_inf)//�����ֱȵ�ǰmax_inf��������ĺ�ѡ�أ����query_result����������ѹ���µĺ�ѡ��
		{
			max_inf = candidate_inf_map[*iter_candidate];
			query_result.clear();
			query_result.push_back(*iter_candidate);
		}
		else if (candidate_inf_map[*iter_candidate] == max_inf)//����ǰ��ѡ����Ӱ����û������뵱ǰ���ֵ��ͬ���򽫷�����ѡ��ѹ��
			query_result.push_back(*iter_candidate);//����ǰ��ѡ��Ӱ���û��������ڵ�ǰ���ֵ����ʲô�������������鿴�µĺ�ѡ��
	}

/*	std::map<geo_coordinate, std::vector<std::string>>::const_iterator iter_result = cand_info.cbegin();
	for (; iter_result != cand_info.cend(); iter_result++)
		std::cout << iter_result->first.get_coordinate().get<0>() << "    " << iter_result->first.get_coordinate().get<1>() << "     " << iter_result->second .size()<< std::endl;
*/

	candidate_inf_map.clear();
	std::vector<geo_coordinate>::const_iterator iter_best = query_result.cbegin();

	std::map<geo_coordinate, std::vector<std::string>>::const_iterator max_coordinate = cand_info.find(*iter_best);
	std::vector<std::string> info = max_coordinate->second;
	candidate_inf_map[*iter_best] = max_coordinate->second.size();
	cand_info.erase(max_coordinate);
	int i = 1;
	while (i < k)
	{
	
		max_inf = 0;
		std::map<geo_coordinate, std::vector<std::string>>::const_iterator iter_cand = cand_info.cbegin();
		for (; iter_cand != cand_info.cend(); iter_cand++)
		{
			std::vector<std::string>::const_iterator iter_info = info.cbegin();
			std::vector<std::string> cand_in = iter_cand->second;
			for (; iter_info != info.cend() && info.size(); iter_info++)
			{
				std::vector<std::string>::const_iterator iter_in = find(cand_in.cbegin(), cand_in.cend(),*iter_info);
				if (iter_in != cand_in.cend() && cand_in.size())
				{
					cand_in.erase(iter_in);//ɾ����info������ͬ���û�
				}
			}
			//std::cout << i << std::endl;
			cand_info[iter_cand->first] = cand_in;//���µ�info��Ϣ�����¸���cand_info
			std::map<geo_coordinate, std::vector<std::string>>::const_iterator iter = cand_info.cbegin();
	/*		for (; iter != cand_info.cend(); iter++)
			{
				std::cout << "see out     " << iter->first.get_coordinate().get<0>() << "     " << iter->first.get_coordinate().get<1>() <<"     "<<iter->second.size()<< std::endl;	
			}*/
		//	iter_cand->second = cand_in;
			if (iter_cand->second.size() > max_inf)
			{
				max_inf = iter_cand->second.size();
				max_coordinate = iter_cand;
			}	
		}
		info = max_coordinate->second;
		candidate_inf_map[max_coordinate->first] = info.size();
		std::cout << max_coordinate->first.get_coordinate().get<0>() << "    " << max_coordinate->first.get_coordinate().get<1>() << "     " << info.size() << std::endl;

		cand_info.erase(max_coordinate);
		i++;
	}
	std::cout << "*****************************************" << std::endl;
	std::map<geo_coordinate, unsigned>::const_iterator iter_map = candidate_inf_map.cbegin();
	for (; iter_map != candidate_inf_map.cend(); iter_map++)
		std::cout << iter_map->first.get_coordinate().get<0>() << "    " << iter_map->first.get_coordinate().get<1>() << "     " << iter_map->second << std::endl;
}

int pino::HashValue(int x)
{
	int hashvalue = (1235 * x + 6583) % MAX_VALUE;
	return hashvalue;
}

int pino::Index(int hashvalue)
{
	int count = 0;
	int flag = 1;
	while (hashvalue)
	{
		if (hashvalue & flag)
		{
			break;
		}
		hashvalue = hashvalue >> 1;
		count++;
	}

	return count;
}

int pino::CountValue(unsigned int bit_vector)
{
	double fai = 0.77351;
	int flag = 1;
	int num, count = 0;
	bit_vector = ~bit_vector;
	while (bit_vector)
	{
		if (bit_vector & flag)
		{
			break;
		}
		bit_vector = bit_vector >> 1;
		count++;
	}
	num = (int)( pow(2,count) / fai);
	return num;
}

void pino::NA_bitVector(std::map<geo_coordinate, int> &query_result, int k)
{
	std::map<geo_coordinate, int> max_result;
	Produce_bitVector(cand_bit_vec, max_result);
	std::map<geo_coordinate, int>::iterator iter_max = max_result.begin();
	std::map<geo_coordinate, int>::iterator iter_cand = cand_bit_vec.begin();
	query_result[iter_max->first] = iter_max->second;
	int temp=iter_max->second;
	cand_bit_vec.erase(iter_max->first);
	for (int i = 1; i < k; i++)
	{
		int max_value = 0;
		for (iter_cand = iter_max = cand_bit_vec.begin(); iter_cand != cand_bit_vec.end(); iter_cand++)
		{
			
			int value = (temp | iter_cand->second) - iter_cand->second;
			if (value >= max_value)
			{
				iter_max = iter_cand;
				max_value = value;
			}
		}
		query_result[iter_max->first] = max_value;
		temp = iter_max->second | temp;
		cand_bit_vec.erase(iter_max->first);
 	}
	
	std::cout << "*****************************************" << std::endl;
	std::map<geo_coordinate, int>::iterator iter_map = query_result.begin();
	for (; iter_map != query_result.cend(); iter_map++)
		std::cout << iter_map->first.get_coordinate().get<0>() << "    " << iter_map->first.get_coordinate().get<1>() << "     " << CountValue(iter_map->second) << std::endl;
	
}

void pino::Produce_bitVector(std::map<geo_coordinate, int> &cand_bit_vec, std::map<geo_coordinate, int> &max_result)
{
	unsigned max_inf = 0;
	int max = 0;
	std::vector<geo_coordinate> query_result;
	// Init inf value for each candidate.
	std::map<geo_coordinate, unsigned> candidate_inf_map;
	std::vector<geo_coordinate>::const_iterator iter_candidate = candidate_vector.cbegin();
	std::map<std::string, std::vector<geo_point>>::const_iterator iter_checkin = checkin_map.cbegin();
	int traj = 0, i = 0;
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate)
	{
		candidate_inf_map[*iter_candidate] = 0;//��ÿ����ѡ��Ӱ����û�������ʼ��Ϊ0
		cand_bit_vec[*iter_candidate] = 0;
	}
	for (; iter_checkin != checkin_map.cend(); iter_checkin++)
	{
		traj_value[iter_checkin->first] = ++i;//Ϊÿһ���û�ָ����ͬ��ֵ
	}

	iter_candidate = candidate_vector.cbegin();//������ѡ��vector
	for (; iter_candidate != candidate_vector.cend(); ++iter_candidate) // Iterate each candidate to calculate its inf value.
	{
//		std::map<std::string, std::vector<geo_point>>::const_iterator iter_checkin = checkin_map.cbegin();

		for (iter_checkin = checkin_map.cbegin(); iter_checkin != checkin_map.cend(); ++iter_checkin) // Iterate each user.
		{
			long double product = 1.0; // "1 - product >= threshold".
			std::vector<geo_point>::const_iterator iter_coordinate = iter_checkin->second.cbegin();
			for (; iter_coordinate != iter_checkin->second.cend(); ++iter_coordinate) // Iterate each position of a user.
			{
				long double distance = geo_distance(iter_candidate->get_coordinate(), *iter_coordinate);
				product *= 1.0 - probability_function(distance);//����Prc��p1����Prc��p2����...��Prc��pn��
			}
			if (1.0 - product >= threshold)//Prc��O��
			{

				++(candidate_inf_map[*iter_candidate]);
				//vec_info.push_back(iter_checkin->first);//��na_k�м����ѡ��Ӱ����û��켣
				int traj = traj_value[iter_checkin->first];
				int hashvalue = HashValue(traj);
				int index = Index(hashvalue);
				cand_bit_vec[*iter_candidate] = cand_bit_vec[*iter_candidate] | (1 << index);
				int a = 1;
			}
		} // Now, all users are checked for this candidate.
		int value = CountValue(cand_bit_vec[*iter_candidate]);
		if (value > max)
		{
			max = value;
			max_result.clear();
			max_result[*iter_candidate] = cand_bit_vec[*iter_candidate];//�ҵ��������㷨��Ӱ���û������Ǹ���ѡ��
		}
		// Update global max inf & current optimal candidate(s).
		if (candidate_inf_map[*iter_candidate] > max_inf)//�����ֱȵ�ǰmax_inf��������ĺ�ѡ�أ����query_result����������ѹ���µĺ�ѡ��
		{
			max_inf = candidate_inf_map[*iter_candidate];
			query_result.clear();
			query_result.push_back(*iter_candidate);
		}
		else if (candidate_inf_map[*iter_candidate] == max_inf)//����ǰ��ѡ����Ӱ����û������뵱ǰ���ֵ��ͬ���򽫷�����ѡ��ѹ��
			query_result.push_back(*iter_candidate);//����ǰ��ѡ��Ӱ���û��������ڵ�ǰ���ֵ����ʲô�������������鿴�µĺ�ѡ��
	}
	int a=2;
}


int pino::pin_vo_k(std::map<geo_coordinate, int> &result, int k)
{
	std::vector<geo_coordinate> query_result;
    query_result.clear();
	int max_inf = 0; // maxminInf.
	cand_info.clear();
	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size(); // Candidates count.

	//auto begin_time = std::chrono::high_resolution_clock::now();
	//auto end_time = std::chrono::high_resolution_clock::now();

	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // The number of candidates inside IA.
		if (iter_user->second.get<U_INF_VALID>()) // Check IA is valid or not.
		{
			// Query candidates inside IA.
			//begin_time = std::chrono::high_resolution_clock::now();
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),
				std::inserter(inf_candidates, inf_candidates.begin()));
			//end_time = std::chrono::high_resolution_clock::now();
			pruned_by_IA += inf_size; // Accumulate inf values.

			// Iterate each candidate in IA.
			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment inf value by 1 for each candidate inside IA.
				cand_info[*iter_inf].push_back(iter_user->first);//���µ�ǰ�û�IA��Χ�ں�ѡ�أ���¼ÿ����ѡ��Ӱ����û��켣
			}//������IA��Χ�ڵĺ�ѡ�أ�ֻ��Ҫ��minInf+1����
		}

		// Query potential candidates inside NIB.
		std::vector<geo_coordinate> potential_candidates;
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		pruned_by_NIB += candidates_count - potential_size; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{
			if (inf_size != 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
				continue;

			candidate_tuple_map[potential_candidates[i]].get<C_TO_VERIFY>().push_back(iter_user->first); // Append this user as a potential user for the candidate.
		}//candidate_tuple_map�����ѡ����Ҫ��֤�û�����Ϣ�����仰˵������NIB��ĺ�ѡ�ؼ���candidate_tuple_map��
	}
	//long long time_cand_counts = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();
	//out_file << "rtree time=" << time_cand_counts << "ns\n";

	// Use a max_heap to rank candidates by potential inf and min inf.
/*	std::priority_queue<CandidateIter> max_heap;
	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		unsigned p_value = iter_tuple->second.get<C_POTENTIAL>() = iter_tuple->second.get<C_INF>() + iter_tuple->second.get<C_TO_VERIFY>().size();	// Potential inf.
		max_heap.push(CandidateIter(p_value, iter_tuple->second.get<C_INF>(), iter_tuple));//�Զѽ��г�ʼ�������Ҹ�����maxInf��ֵ
	}

	//out_file << "New PIN-VO\n";
	//begin_time = std::chrono::high_resolution_clock::now();
	// Iterate candidates based on (potential and min infs) max heap.����Ҫ��֤�ĺ�ѡ�ؽ��д���
	while (!max_heap.empty())
	{

		geo_point candidate = max_heap.top().iter_map->first.get_coordinate(); // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = max_heap.top().iter_map->second.get<C_TO_VERIFY>().cbegin(),
			iter_end = max_heap.top().iter_map->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
				iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate, *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(max_heap.top().iter_map->second.get<C_INF>()); // Update "minInf" value.
					cand_info[candidate].push_back(*iter_user);
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}
		}
		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);

		max_heap.pop(); // Pop next candidate to heap top.
	}
	*/


	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		geo_coordinate candidate = iter_tuple->first; // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = iter_tuple->second.get<C_TO_VERIFY>().cbegin(),
			iter_end = iter_tuple->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
				iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate.get_coordinate(), *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(iter_tuple->second.get<C_INF>()); // Update "minInf" value.
					cand_info[candidate].push_back(*iter_user);
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}
		}
		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);
		
	}



	int count = 0;

	std::map<geo_coordinate, std::vector<std::string>>::iterator  iter_cand,iter_max;
//	
	std::vector<std::string> info,cand_in;
	std::vector<std::string>::iterator iter_info;
	std::vector<geo_coordinate>::iterator iter_result = query_result.begin();
	info = cand_info[*iter_result];
	result[*iter_result] = max_inf;
	cand_info.erase(*iter_result);
	for (int i = 1; i < k; i++)
	{ 
		max_inf = 0;

		for (iter_cand = cand_info.begin(); iter_cand != cand_info.end(); iter_cand++)
		{
			cand_in = iter_cand->second;
			for (iter_info = info.begin(); iter_info != info.end() && info.size(); iter_info++)
			{
				std::vector<std::string>::const_iterator iter_in = find(cand_in.cbegin(), cand_in.cend(), *iter_info);
				if (iter_in != cand_in.cend() && cand_in.size())
				{
					cand_in.erase(iter_in);//ɾ����info������ͬ���û�
				}
			}
			cand_info[iter_cand->first] = cand_in;
			if (cand_info[iter_cand->first].size() > max_inf)
			{
				max_inf = cand_info[iter_cand->first].size();
				iter_max = iter_cand;
			}
		}
		info = iter_max->second;
		result[iter_max->first] = max_inf;
	//	std::cout << iter_max->first.get_coordinate().get<0>() << "    " << iter_max->first.get_coordinate().get<1>() << "     " << max_inf << std::endl;

		cand_info.erase(iter_max->first);
	}
	std::cout << "****************************************" << std::endl;
	std::map<geo_coordinate, int>::iterator iter;
	for (iter = result.begin(); iter != result.end(); iter++)
	{
		count += iter->second;
		std::cout << iter->first.get_coordinate().get<0>() << "    " << iter->first.get_coordinate().get<1>() << "     " << iter->second << std::endl;
	}
	std::cout << count << std::endl; 
	return count;
}

void pino::Pin_Vo_vitVector(std::map<geo_coordinate, int> &result, int k)
{
	std::vector<geo_coordinate> query_result;
	query_result.clear();
	int max_inf = 0, max = 0 ; // maxminInf.
	cand_info.clear();
	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size(); // Candidates count.
	std::map < std::string, std::vector<geo_point>>::const_iterator iter_checkin;
	//auto begin_time = std::chrono::high_resolution_clock::now();
	//auto end_time = std::chrono::high_resolution_clock::now();
	int i=0;
	for (iter_checkin = checkin_map.cbegin(); iter_checkin != checkin_map.cend(); iter_checkin++)
	{
		traj_value[iter_checkin->first] = ++i;//Ϊÿһ���û�ָ����ͬ��ֵ
	}
	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // The number of candidates inside IA.
		if (iter_user->second.get<U_INF_VALID>()) // Check IA is valid or not.
		{
			// Query candidates inside IA.
			//begin_time = std::chrono::high_resolution_clock::now();
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),
				std::inserter(inf_candidates, inf_candidates.begin()));
			//end_time = std::chrono::high_resolution_clock::now();
			pruned_by_IA += inf_size; // Accumulate inf values.

			// Iterate each candidate in IA.
			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment inf value by 1 for each candidate inside IA.
				//cand_info[*iter_inf].push_back(iter_user->first);//���µ�ǰ�û�IA��Χ�ں�ѡ�أ���¼ÿ����ѡ��Ӱ����û��켣
				int trj = traj_value[iter_user->first];
				int hashvalue = HashValue(trj);
				int index = Index(hashvalue);
				
				cand_bit_vec[*iter_inf] = cand_bit_vec[*iter_inf] | (1 << index);
			}//������IA��Χ�ڵĺ�ѡ�أ�ֻ��Ҫ��minInf+1����
		}

		// Query potential candidates inside NIB.
		std::vector<geo_coordinate> potential_candidates;
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		pruned_by_NIB += candidates_count - potential_size; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{
			if (inf_size != 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
				continue;

			candidate_tuple_map[potential_candidates[i]].get<C_TO_VERIFY>().push_back(iter_user->first); // Append this user as a potential user for the candidate.
		}//candidate_tuple_map�����ѡ����Ҫ��֤�û�����Ϣ�����仰˵������NIB��ĺ�ѡ�ؼ���candidate_tuple_map��
	}
	//long long time_cand_counts = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();
	//out_file << "rtree time=" << time_cand_counts << "ns\n";

	// Use a max_heap to rank candidates by potential inf and min inf.
/*	std::priority_queue<CandidateIter> max_heap;
	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		unsigned p_value = iter_tuple->second.get<C_POTENTIAL>() = iter_tuple->second.get<C_INF>() + iter_tuple->second.get<C_TO_VERIFY>().size();	// Potential inf.
		max_heap.push(CandidateIter(p_value, iter_tuple->second.get<C_INF>(), iter_tuple));//�Զѽ��г�ʼ�������Ҹ�����maxInf��ֵ
	}

	//out_file << "New PIN-VO\n";
	//begin_time = std::chrono::high_resolution_clock::now();
	// Iterate candidates based on (potential and min infs) max heap.����Ҫ��֤�ĺ�ѡ�ؽ��д���
	while (!max_heap.empty())
	{
		// Upper and lower bounds.
		if (max_heap.top().potential_value < max_inf)
			break;
		//out_file << max_heap.top().iter_map->first.longitude() << ", " << max_heap.top().iter_map->first.latitude() << ", g=" << max_inf << ", p=" << max_heap.top().potential_value << ", i=" << max_heap.top().min_inf_value << "\n";

		geo_point candidate = max_heap.top().iter_map->first.get_coordinate(); // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = max_heap.top().iter_map->second.get<C_TO_VERIFY>().cbegin(),
			iter_end = max_heap.top().iter_map->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
				iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate, *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(max_heap.top().iter_map->second.get<C_INF>()); // Update "minInf" value.
			//		cand_info[candidate].push_back(*iter_user);
					int trj = traj_value[*iter_user];
					int hashvalue = HashValue(trj);
					int index = Index(hashvalue);
					cand_bit_vec[candidate] = cand_bit_vec[candidate] | ( 1 << index);
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}

			if (!is_inf // Not inf this user.
				&& --(max_heap.top().iter_map->second.get<C_POTENTIAL>()) < max_inf) // Potential <= global max inf.
				break; // The candidate is not optimal.�����ѡ�ز���Ӱ����û�����maxInf-1�����������1����ֹͣ�Ըú�ѡ�ص���֤
		}

		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);

		max_heap.pop(); // Pop next candidate to heap top.
	}*/

	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		geo_coordinate candidate = iter_tuple->first; // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = iter_tuple->second.get<C_TO_VERIFY>().cbegin(),
			iter_end = iter_tuple->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
				iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate.get_coordinate(), *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(iter_tuple->second.get<C_INF>()); // Update "minInf" value.
				//	cand_info[candidate].push_back(*iter_user);
					int trj = traj_value[*iter_user];
					int hashvalue = HashValue(trj);
					int index = Index(hashvalue);
					cand_bit_vec[candidate] = cand_bit_vec[candidate] | (1 << index);
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}
		}
		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);

	}



	int current = 0;
	std::map<geo_coordinate, int>::iterator iter_cur = cand_bit_vec.begin(), iter_max;
	for (int i = 0; i < k; i++)
	{
		int max = 0;
	//	iter_max = iter_cur;
		for (iter_cur = iter_max = cand_bit_vec.begin(); iter_cur != cand_bit_vec.end(); iter_cur++)
		{
			int value;
			if (i == 0)
				value = current | iter_cur->second;
			else
				value = (current | iter_cur->second) - current;
		//	iter_cur->second = value;
			if (CountValue(value) >= CountValue(max))
			{
				max = value;
				iter_max = iter_cur;
			}
			int a = 3;
		}
	/*	if (i == 0)
			current = max;*/
		current = current | max;
		result[iter_max->first] = max;
		std::cout << iter_max->first.get_coordinate().get<0>() << "    " << iter_max->first.get_coordinate().get<1>() << "     " << CountValue(max) << std::endl;

		cand_bit_vec.erase(iter_max->first);
		int a = 1;
	} 
	std::cout << "****************************************" << std::endl;
	std::map<geo_coordinate, int>::iterator iter;
	for (iter = result.begin(); iter != result.end(); iter++)
	{
		std::cout << iter->first.get_coordinate().get<0>() << "    " << iter->first.get_coordinate().get<1>() << "     " << CountValue(iter->second) << std::endl;
	}
}

void pino::init()
{
	srand((unsigned)time(NULL));
	candidate_inf.clear();
	std::cout << "parameter:" << std::endl;
	for (int i = 0; i < HASH_NUM; i++)
	{
		function[i] = random(1000,5000);
		
	//	std::cout << function[i] << ", " << std::endl;
	}
	std::vector<geo_coordinate>::iterator it = candidate_vector.begin();
	for (it; it != candidate_vector.end(); it++)
	{
		for (int i = 0; i < HASH_NUM; i++)
			candidate_inf[*it].push_back(0);
		
	}
	std::map < std::string, std::vector<geo_point>>::const_iterator iter_checkin;
	int i = 0;
	for (iter_checkin = checkin_map.cbegin(); iter_checkin != checkin_map.cend(); iter_checkin++)
	{
		traj_value[iter_checkin->first] = ++i;//Ϊÿһ���û�ָ����ͬ��ֵ
	}

}
void pino::new_hash(geo_coordinate point, int value)
{
	for (int i = 0; i < HASH_NUM; i++)
	{
		int hash_result = (function[i] * value + 2361) % MAX_VALUE;
		int hash_index = Index(hash_result);
		candidate_inf[point][i] = candidate_inf[point][i] | (1 << hash_index);
	}
	int a = 0;
}
int pino::num_one(unsigned int value)
{
	int count = 0;
	while (value)
	{
/*		if (value & 1)
		{
			count++;
			value = value >> 1;
		}
		else
			break;*/
		if (value)
		{
			if (value & 1)
				count++;
			value = value >> 1;
		}

	}
	return count;
}
bool cmp(pino::Node &x, pino::Node &y)
{
	if (x.count> y.count)
		return true;
	else
		return false;
}
int pino::new_index(unsigned value)
{
	return value % 32;
}
void pino::pin_vo_num_vitVector(std::map<geo_coordinate, int> &result, int k)
{
	init();
	std::vector<geo_coordinate> query_result;
	query_result.clear();
	int max_inf = 0; // maxminInf.
	cand_info.clear();
	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size(); // Candidates count.

	//auto begin_time = std::chrono::high_resolution_clock::now();
	//auto end_time = std::chrono::high_resolution_clock::now();


	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // The number of candidates inside IA.
		if (iter_user->second.get<U_INF_VALID>()) // Check IA is valid or not.
		{
			// Query candidates inside IA.
			//begin_time = std::chrono::high_resolution_clock::now();
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),
				std::inserter(inf_candidates, inf_candidates.begin()));
			//end_time = std::chrono::high_resolution_clock::now();
			pruned_by_IA += inf_size; // Accumulate inf values.

			// Iterate each candidate in IA.
			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment inf value by 1 for each candidate inside IA.
			//	cand_info[*iter_inf].push_back(iter_user->first);//���µ�ǰ�û�IA��Χ�ں�ѡ�أ���¼ÿ����ѡ��Ӱ����û��켣
				int traj = traj_value[iter_user->first];
				new_hash(*iter_inf,traj);//���ú�������ù켣�Ĳ�ͬ��ϣ����ֵ
			}//������IA��Χ�ڵĺ�ѡ�أ�ֻ��Ҫ��minInf+1����
		}

		// Query potential candidates inside NIB.
		std::vector<geo_coordinate> potential_candidates;
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		pruned_by_NIB += candidates_count - potential_size; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{
			if (inf_size != 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
				continue;

			candidate_tuple_map[potential_candidates[i]].get<C_TO_VERIFY>().push_back(iter_user->first); // Append this user as a potential user for the candidate.
		}//candidate_tuple_map�����ѡ����Ҫ��֤�û�����Ϣ�����仰˵������NIB��ĺ�ѡ�ؼ���candidate_tuple_map��
	}


	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		geo_coordinate candidate = iter_tuple->first; // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = iter_tuple->second.get<C_TO_VERIFY>().cbegin(),
			iter_end = iter_tuple->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
				iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate.get_coordinate(), *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(iter_tuple->second.get<C_INF>()); // Update "minInf" value.
	//				cand_info[candidate].push_back(*iter_user);
					int traj = traj_value[*iter_user];
					new_hash(candidate,traj);
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}
		}
		if (min_inf > max_inf) // Update global (max_inf) maxminInf.
		{
			max_inf = min_inf;
			query_result.clear();
			query_result.push_back(candidate);
		}
		else if (min_inf == max_inf)
			query_result.push_back(candidate);

	}

//	std::priority_queue<Node> max_heap;
	std::vector<Node> candidate_bit;
	std::map<geo_coordinate, std::vector<unsigned int>>::iterator iter_cand_inf = candidate_inf.begin();
	for (iter_cand_inf; iter_cand_inf != candidate_inf.end(); iter_cand_inf++)
	{
		unsigned int sum = 0;
		for (int i = 0; i < HASH_NUM; i++)
			sum =sum + num_one(candidate_inf[iter_cand_inf->first][i]);
	//	candidate_bit.pushback(Node(iter_cand_inf->first,sum,&iter_cand_inf->second));
		candidate_bit.push_back(Node(iter_cand_inf->first, sum, iter_cand_inf->second));
	}
	//std::sort(candidate_bit.begin(),candidate_bit.end(),cmp);
	std::sort(candidate_bit.begin(),  candidate_bit.end(),cmp);
	int a = 0;
	std::vector<unsigned int> current,temp;
	for (int i = 0; i < HASH_NUM; i++)
		current.push_back(0);
	for (int i = 0; i < k ; i++)
	{
		unsigned int max = 0;
		std::vector<Node>::iterator iter_cur,iter_max;
		for (iter_cur = iter_max = candidate_bit.begin(); iter_cur != candidate_bit.end(); iter_cur++)
		{
			if (max > iter_cur->count)
				break;
			unsigned int value = 0;
			for (int j = 0; j < HASH_NUM; j++)
			{
				value += num_one(((*iter_cur).bit_map[j] | current[j]) - current[j]);
			}
		//	std::cout << value <<"-------------"<< max << std::endl;
			if (max < value)
			{
				max = value;
				iter_max = iter_cur;
			}
		}
		for (int j = 0; j < HASH_NUM; j++)
		{
			current[j] = current[j] | (*iter_max).bit_map[j];
			int a = 1;
		}
		result[iter_max->candidate] = max;
		candidate_bit.erase(iter_max);	
	//	std::cout << i << std::endl;
	}
	std::cout << "****************************************" << std::endl;
	std::map<geo_coordinate, int>::iterator iter;
	for (iter = result.begin(); iter != result.end(); iter++)
	{
	//	std::cout << iter->first.get_coordinate().get<0>() << "    " << iter->first.get_coordinate().get<1>() << "     " << iter->second << std::endl;
		std::cout << iter->first.get_coordinate().get<0>() << "    " << iter->first.get_coordinate().get<1>() << std::endl;

	}
}

int pino::test_traj(std::map<geo_coordinate, int> result,int k)
{
	int count = 0;
	cand_info.clear();
	std::map<geo_coordinate, int>::iterator iter_result;
	std::map<std::string, std::vector<geo_point>>::const_iterator iter_user;
	for (iter_result = result.begin(); iter_result != result.end(); iter_result++)
	{
		for (iter_user = checkin_map.cbegin(); iter_user != checkin_map.cend(); iter_user++)
		{
			long double product = 1.0; // "1 - product >= threshold".
			std::vector<geo_point>::const_iterator iter_coordinate = iter_user->second.cbegin();
			for (; iter_coordinate != iter_user->second.cend(); ++iter_coordinate) // Iterate each position of a user.
			{
				long double distance = geo_distance(iter_result->first.get_coordinate(), *iter_coordinate);
				product *= 1.0 - probability_function(distance);//����Prc��p1����Prc��p2����...��Prc��pn��
			}
			if (1.0 - product >= threshold)//Prc��O��
				cand_info[iter_result->first].push_back(iter_user->first);
		}		
	}
	std::map<geo_coordinate, std::vector<std::string>>::iterator iter_cand, iter_max;
	std::vector<std::string> info, cand_in;
	std::vector<std::string>::iterator iter_info;
	for (int i = 0; i < k; i++)
	{
		int max = -1;
		if (i)
		{
			
			for (iter_cand = cand_info.begin(); iter_cand != cand_info.end(); iter_cand++)
			{
				cand_in = iter_cand->second;
				for (iter_info = info.begin(); iter_info != info.end(); iter_info++)
				{
					std::vector<std::string>::const_iterator iter_in = find(cand_in.cbegin(), cand_in.cend(), *iter_info);
					if (iter_in != cand_in.cend() && cand_in.size())
					{
						cand_in.erase(iter_in);//ɾ����info������ͬ���û�
					}
				}
				cand_info[iter_cand->first] = cand_in;
			}
		}
		for (iter_cand = iter_max =cand_info.begin(); iter_cand != cand_info.end(); iter_cand++)
		{
	//		std::cout << iter_cand->second.size() << std::endl;
			if (iter_cand->second.size() > max)
			{
				max = iter_cand->second.size();
				iter_max = iter_cand;
			}
		}
		info = iter_max->second;
		std::cout << iter_max->first.get_coordinate().get<0>() << "    " << iter_max->first.get_coordinate().get<1>() << "    " << info.size() << std::endl;
		count += info.size();
		cand_info.erase(iter_max);
	}
	return count;
}

void pino::Baseline(std::vector<geo_coordinate> &result, int k)
{
//	std::vector<geo_coordinate> query_result;
//	query_result.clear();
	cand_info.clear();
	std::map<int, geo_coordinate> Order_cand;
//	int max_inf = 0; // maxminInf.
	cand_info.clear();
	pruned_by_IA = pruned_by_NIB = 0; // Init the two number.
	unsigned candidates_count = candidate_tuple_map.size(); // Candidates count.

	std::map<std::string, user_box_tuple>::const_iterator iter_user = user_box_tuple_map.cbegin();
	for (; iter_user != user_box_tuple_map.cend(); ++iter_user) // Iterate each user.
	{
		std::set<geo_coordinate> inf_candidates; // Candidates inside IA.
		unsigned inf_size = 0; // The number of candidates inside IA.
		if (iter_user->second.get<U_INF_VALID>()) // Check IA is valid or not.
		{
			// Query candidates inside IA.
			//begin_time = std::chrono::high_resolution_clock::now();
			inf_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_INF_BOUND>()),
				std::inserter(inf_candidates, inf_candidates.begin()));
			//end_time = std::chrono::high_resolution_clock::now();
			pruned_by_IA += inf_size; // Accumulate inf values.

			// Iterate each candidate in IA.
			for (std::set<geo_coordinate>::const_iterator iter_inf = inf_candidates.cbegin(); iter_inf != inf_candidates.cend(); ++iter_inf)
			{
				++(candidate_tuple_map[*iter_inf].get<C_INF>()); // Increment inf value by 1 for each candidate inside IA.
				cand_info[*iter_inf].push_back(iter_user->first);//���µ�ǰ�û�IA��Χ�ں�ѡ�أ���¼ÿ����ѡ��Ӱ����û��켣
			}//������IA��Χ�ڵĺ�ѡ�أ�ֻ��Ҫ��minInf+1����
		}

		// Query potential candidates inside NIB.
		std::vector<geo_coordinate> potential_candidates;
		unsigned potential_size = candidate_rtree.query(boost::geometry::index::intersects(iter_user->second.get<U_POTENTIAL>()), std::back_inserter(potential_candidates));
		pruned_by_NIB += candidates_count - potential_size; // Accumulate non-inf values.

		for (unsigned i = 0; i < potential_size; ++i) // Iterate each potential candidate.
		{
			if (inf_size != 0 // IA and IA candidates exist.
				&& inf_candidates.find(potential_candidates[i]) != inf_candidates.end()) // The candidate is exactly on the IA MBR bound. No need to check it.
				continue;

			candidate_tuple_map[potential_candidates[i]].get<C_TO_VERIFY>().push_back(iter_user->first); // Append this user as a potential user for the candidate.
		}//candidate_tuple_map�����ѡ����Ҫ��֤�û�����Ϣ�����仰˵������NIB��ĺ�ѡ�ؼ���candidate_tuple_map��
	}

	std::map<geo_coordinate, candidate_tuple>::iterator iter_tuple = candidate_tuple_map.begin();
	for (; iter_tuple != candidate_tuple_map.end(); ++iter_tuple)
	{
		geo_coordinate candidate = iter_tuple->first; // (Heap top) Candidate coordinate.
		unsigned min_inf = 0; // Actual "minInf" value.

		std::vector<std::string>::const_iterator iter_user = iter_tuple->second.get<C_TO_VERIFY>().cbegin(),
			iter_end = iter_tuple->second.get<C_TO_VERIFY>().cend();
		//���ò��Զ�
		for (; iter_user != iter_end; ++iter_user) // Iterate potential users' strings.
		{
			std::vector<geo_point>::const_iterator iter_pos = checkin_map[*iter_user].cbegin(),
				iter_pos_end = checkin_map[*iter_user].cend();
			long double product = 1.0;
			bool is_inf = false;
			for (; iter_pos != iter_pos_end; ++iter_pos) // Iterate positions of a potential user.
			{
				long double distance = geo_distance(candidate.get_coordinate(), *iter_pos);
				product *= 1.0 - probability_function(distance);

				if (1.0 - product >= threshold) // Early stopping.
				{
					is_inf = true;
					min_inf = ++(iter_tuple->second.get<C_INF>()); // Update "minInf" value.
					cand_info[candidate].push_back(*iter_user);
					break; // Jump out and compute next user (next iter_user).������ѭ����ת��ȥ��֤��һ���û�
				}
			}
		}
	}
	std::map<geo_coordinate, std::vector<std::string>>::iterator iter_cand = cand_info.begin();
	for (iter_cand; iter_cand != cand_info.end(); iter_cand++)
		Order_cand[iter_cand->second.size()] = iter_cand->first;
	std::map<int, geo_coordinate>::iterator iter_order = Order_cand.end();
	iter_order--;
	for (int i = 0; i < k; i++)
	{
		result.push_back(iter_order->second);
		iter_order--;
	}
}

int pino::test_traj1(std::vector<geo_coordinate> result,int k)
{
	int count = 0;
//	cand_info.clear();
	std::vector<geo_coordinate>::iterator iter_result;
	std::map<std::string, std::vector<geo_point>>::const_iterator iter_user;
	std::map<geo_coordinate, std::vector<std::string>> new_cand;
	std::map<geo_coordinate, std::vector<std::string>>::iterator iter_cand, iter_max;
	std::vector<std::string> info, cand_in;
	std::vector<std::string>::iterator iter_info;
	for (iter_result = result.begin(); iter_result != result.end(); iter_result++)
	{
		new_cand[*iter_result] = cand_info[*iter_result];
	}
	for (int i = 0; i < k; i++)
	{
		int max = -1;
		if (i)
		{

			for (iter_cand = new_cand.begin(); iter_cand != new_cand.end(); iter_cand++)
			{
				cand_in = iter_cand->second;
				for (iter_info = info.begin(); iter_info != info.end(); iter_info++)
				{
					std::vector<std::string>::const_iterator iter_in = find(cand_in.cbegin(), cand_in.cend(), *iter_info);
					if (iter_in != cand_in.cend() && cand_in.size())
					{
						cand_in.erase(iter_in);//ɾ����info������ͬ���û�
					}
				}
				new_cand[iter_cand->first] = cand_in;
			}
		}
		for (iter_cand = iter_max = new_cand.begin(); iter_cand != new_cand.end(); iter_cand++)
		{
			//		std::cout << iter_cand->second.size() << std::endl;
			if (iter_cand->second.size() > max)
			{
				max = iter_cand->second.size();
				iter_max = iter_cand;
			}
		}
		info = iter_max->second;
		std::cout << iter_max->first.get_coordinate().get<0>() << "    " << iter_max->first.get_coordinate().get<1>() << "    " << info.size() << std::endl;
		count += info.size();
		new_cand.erase(iter_max);
	}
	return count;
	return count;
}
