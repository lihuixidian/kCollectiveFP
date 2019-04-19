#pragma once

#include "stdafx.h"

typedef boost::geometry::model::point<long double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree>> geo_point;
typedef boost::geometry::model::box<geo_point> geo_box;//是一个矩形
typedef boost::geometry::model::multi_point<geo_point> geo_multi_point;
typedef boost::geometry::model::polygon<geo_point> geo_polygon;//多边形
typedef boost::geometry::model::ring<geo_point> geo_ring;
/*这里的ring就是我们通常说的多边形闭合区域(内部不存在缕空)，模板参数为true，表示顺时针存储点，为false，
表示逆时针存储点，由于MM_TEXT坐标系与传统上的坐标系的Y轴方向是相反的，所以最后为false，
将TopLeft、TopRight、BottomRight、BottomLeft、TopLeft以此存储到ring中，以便能正确计算  */
typedef boost::geometry::model::segment<geo_point> geo_segment;//线段
typedef boost::geometry::model::linestring<geo_point> geo_linestring;//由多个点构成的组合，也可以表示多边形
typedef boost::geometry::index::rtree<geo_point, boost::geometry::index::rstar<8>> geo_point_rtree;