#pragma once

#include <string>
#include <ctime>

#include "typedefs.h"

#include "geo_coordinate.h"


class checkin//记录用户每次登录的地点，时间
{
public:
	enum checkin_time_format { citf_all,		// All fields.
							   citf_hour_min };	// hour and minute fields.登记时间格式

public:
	checkin(){};
	//若是同一个用户，比较年份。若年份相同，则比较月份。。。以此类推，<则返回true
	checkin(const std::string& user, long double longitude, long double latidute, checkin_time_format time_format, const std::tm& time);
	virtual ~checkin(){};

public:
	bool operator<(const checkin& rhs) const;

// Attributes
public:
	long double longitude() const { return checkin_coordinate.longitude(); };
	long double latitude() const { return checkin_coordinate.latitude(); };
	void longitude(long double longitude) { checkin_coordinate.longitude(longitude); };
	void latitude(long double latitude) { checkin_coordinate.latitude(latitude); };

// Members
public:
	std::string user;
	geo_coordinate checkin_coordinate;
	checkin_time_format time_format;
	std::tm time; 
};