#pragma once
#include<bits/stdc++.h>
//邻接表
#define MAX_NEIGHBOR_COUNTS 80
class NeighborTable {
public:
	NeighborTable();
	void reset(unsigned short pointCounts);//重置邻接表
	void point_prepare(unsigned short ptIndex);//预备点数据
	bool point_add_neighbor(unsigned short ptIndex, float distance);//给当前点添加邻接表
	void point_commit();//递交点给邻接表
	int getNeighborCounts(unsigned short ptIndex) { return m_pointExtraData[ptIndex].neighborCounts; }//获取邻接表中的点个数
	void getNeighborInfo(unsigned short ptIndex, int index, unsigned int& neighborIndex, float& neighborDistance);//获取索引ptIndex的邻接表中第index个点的数据

	~NeighborTable();
private:
	union PointExtraData {
		struct {
			unsigned neighborDataOffset : 24;//偏移
			unsigned neighborCounts : 8;//个数
		};
	};
	PointExtraData* m_pointExtraData;//邻接表信息
	unsigned int m_pointCounts;//粒子数
	unsigned int m_pointCapacity;//粒子容量
	unsigned char* m_neighborDataBuf;// 邻接表的数据缓存
	unsigned int m_dataBufSize;//bytes 数据缓存尺寸
	unsigned int m_dataBufOffset;//bytes 数据缓存中的偏移

	//当前点点数据
	unsigned short m_currPoint;//索引
	int m_currNeighborCounts;//邻居点数量
	unsigned short m_currNeightborIndex[MAX_NEIGHBOR_COUNTS];//邻居中点的索引
	float m_currNrighborDistance[MAX_NEIGHBOR_COUNTS];//邻居中点的距离

	void _growDataBuf(unsigned int need_size);//扩容
};

NeighborTable::NeighborTable() :
	m_pointExtraData(0),
	m_pointCounts(0),
	m_pointCapacity(0),
	m_neighborDataBuf(0),
	m_dataBufSize(0),
	m_dataBufOffset(0),
	m_currPoint(0),
	m_currNeighborCounts(0) {
}

NeighborTable::~NeighborTable() {
	if (m_pointExtraData) {
		delete[] m_pointExtraData;
		m_pointExtraData = nullptr;
	}
	if (m_neighborDataBuf) {
		delete[] m_neighborDataBuf;
		m_neighborDataBuf = nullptr;
	}
}

void NeighborTable::reset(unsigned short pointCounts) {
	int a = sizeof(PointExtraData);
	if (pointCounts > m_pointCapacity) {
		if (m_pointExtraData) {
			delete[] m_pointExtraData;
			m_pointExtraData = nullptr;
		}
		m_pointExtraData = new PointExtraData[a * pointCounts]();
		m_pointCapacity = m_pointCounts;
	}
	m_pointCounts = pointCounts;
	memset(m_pointExtraData, 0, a * m_pointCapacity);
	m_dataBufOffset = 0;
}
//准备点
void NeighborTable::point_prepare(unsigned short ptIndex) {
	m_currPoint = ptIndex;
	m_currNeighborCounts = 0;
}

bool NeighborTable::point_add_neighbor(unsigned short ptIndex, float distance) {
	if (m_currNeighborCounts >= MAX_NEIGHBOR_COUNTS)return false;
	m_currNeightborIndex[m_currNeighborCounts] = ptIndex;
	m_currNrighborDistance[m_currNeighborCounts] = distance;
	m_currNeighborCounts++;
	return true;
}
void NeighborTable::point_commit() {
	if (m_currNeighborCounts == 0)return;
	unsigned int index_size = m_currNeighborCounts * sizeof(unsigned short);
	unsigned int distance_size = m_currNeighborCounts * sizeof(float);

	//扩大空间
	if (m_dataBufOffset + index_size + distance_size > m_dataBufSize)
		_growDataBuf(m_dataBufOffset + index_size + distance_size);

	//设置邻居数据
	m_pointExtraData[m_currPoint].neighborCounts = m_currNeighborCounts;
	m_pointExtraData[m_currPoint].neighborDataOffset = m_dataBufOffset;

	//复制索引点的信息到数据缓存中
	memcpy(m_neighborDataBuf + m_dataBufOffset, m_currNeightborIndex, index_size);
	m_dataBufOffset += index_size;
	memcpy(m_neighborDataBuf + m_dataBufOffset, m_currNrighborDistance, distance_size);
	m_dataBufOffset += distance_size;
}

void NeighborTable::_growDataBuf(unsigned int need_size) {
	unsigned int newSize = m_dataBufSize > 0 ? m_dataBufSize : 1;
	while (newSize < need_size) {
		newSize *= 2;
	}
	if (newSize < 2024) {
		newSize = 1024;
	}
	unsigned char* newBuf = new unsigned char[newSize]();
	if (m_neighborDataBuf) {
		memcpy(newBuf, m_neighborDataBuf, m_dataBufSize);
		delete[] m_neighborDataBuf;
	}
	m_neighborDataBuf = newBuf;
	m_dataBufSize = newSize;
}

//ptIndex是点的索引，index是该点邻居表内的邻居索引,获取的neighborIndex为该邻居的默认索引
void NeighborTable::getNeighborInfo(unsigned short ptIndex, int index, unsigned int& neighborIndex, float& neighborDistance) {
	PointExtraData neighData = m_pointExtraData[ptIndex];
	unsigned short* indexBuf = (unsigned short*)(m_neighborDataBuf + neighData.neighborDataOffset);
	float* distanceBuf = (float*)(m_neighborDataBuf + neighData.neighborDataOffset + sizeof(unsigned short) * neighData.neighborCounts);

	neighborIndex = indexBuf[index];
	neighborDistance = distanceBuf[index];
}