#include "stdafx.h"
#include "query_generator.h"


void query_generator::random(const std::vector<geo_coordinate>& src, std::vector<geo_coordinate>& obj, unsigned count)
{
	obj.clear();
	std::vector<geo_coordinate> residual(src);//创建residual容器的时候，传入了src的内容，相当于两个容器中的内容是一样的
	srand(static_cast<unsigned>(time(NULL)));//是初始化随机函数种子
	while (obj.size() < count)
	{
		int random_index = static_cast<int>((double)rand() / (RAND_MAX + 1) * (residual.size() - 1));
		if (find(obj.begin(), obj.end(), residual[random_index]) == obj.end())	// Not the same one.每次找obj数组中没有的元素
			obj.push_back(residual[random_index]);

		// Insert a different coordinate to candidate vector.
		std::vector<geo_coordinate>::iterator iter = find(residual.begin(), residual.end(), residual[random_index]);
		residual.erase(iter);
	}
}
