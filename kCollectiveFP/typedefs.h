#pragma once

#include "stdafx.h"

typedef boost::geometry::model::point<long double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree>> geo_point;
typedef boost::geometry::model::box<geo_point> geo_box;//��һ������
typedef boost::geometry::model::multi_point<geo_point> geo_multi_point;
typedef boost::geometry::model::polygon<geo_point> geo_polygon;//�����
typedef boost::geometry::model::ring<geo_point> geo_ring;
/*�����ring��������ͨ��˵�Ķ���αպ�����(�ڲ��������ƿ�)��ģ�����Ϊtrue����ʾ˳ʱ��洢�㣬Ϊfalse��
��ʾ��ʱ��洢�㣬����MM_TEXT����ϵ�봫ͳ�ϵ�����ϵ��Y�᷽�����෴�ģ��������Ϊfalse��
��TopLeft��TopRight��BottomRight��BottomLeft��TopLeft�Դ˴洢��ring�У��Ա�����ȷ����  */
typedef boost::geometry::model::segment<geo_point> geo_segment;//�߶�
typedef boost::geometry::model::linestring<geo_point> geo_linestring;//�ɶ���㹹�ɵ���ϣ�Ҳ���Ա�ʾ�����
typedef boost::geometry::index::rtree<geo_point, boost::geometry::index::rstar<8>> geo_point_rtree;