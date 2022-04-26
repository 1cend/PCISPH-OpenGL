#pragma once
#include"point.hpp"
//空间网格类
class fBox3 {
public:
	fBox3() :min(glm::vec3(0)), max(glm::vec3(0)) {}
	fBox3(glm::vec3 aMin, glm::vec3 aMax) :min(aMin), max(aMax) {}
	glm::vec3 min, max;
};

class GridContainer {
public:
	GridContainer() {}
	~GridContainer() {}
	//空间细分
	void init(const fBox3& box, float sim_scale, float cell_size, float border, int* rex);
	void insertParticles(PointBuffer* pointBuffer);
	void findCells(const glm::vec3& p, float radius, int* gridCell);
	void findTwoCells(const glm::vec3& p, float radius, int* gridCell);
	int findCell(const glm::vec3& p);
	int getGridData(int gridIndex);

	const glm::ivec3* getGridRes()const { return &m_GridRes; }
	const glm::vec3* getGridMin(void) const { return &m_GridMin; }
	const glm::vec3* getGridMax(void) const { return &m_GridMax; }
	const glm::vec3* getGridSize(void) const { return &m_GridSize; }

	int getGridCellIndex(float px, float py, float pz);//获取对应xyz下的网格索引
	float getDelta() { return m_GridDelta.x; }

private:
	//空间网格
	std::vector<int> m_gridData;//网格信息（储存网格内的当前的粒子）
	glm::vec3 m_GridMin;//网格左下角
	glm::vec3 m_GridMax;//网格右上角
	glm::ivec3 m_GridRes;//网格规格（n * m * l）
	glm::vec3 m_GridSize;//网格大小
	glm::vec3 m_GridDelta;//网格偏移量
	float m_GridCellSize;//一个格子大小（通常为2倍的光滑核半径）
};

int GridContainer::getGridData(int gridIndex) {
	if (gridIndex < 0 || gridIndex >= m_gridData.size())return -1;
	return m_gridData[gridIndex];
}

int GridContainer::getGridCellIndex(float px, float py, float pz) {
	int gx = (int)((px - m_GridMin.x) * m_GridDelta.x);
	int gy = (int)((py - m_GridMin.y) * m_GridDelta.y);
	int gz = (int)((pz - m_GridMin.z) * m_GridDelta.z);
	return (gz * m_GridRes.y + gy) * m_GridRes.x + gx;
}

void GridContainer::init(const fBox3& box, float sim_scale, float cell_size, float border, int* rex) {
	float world_cellsize = cell_size / sim_scale;

	m_GridMin = box.min;
	m_GridMin -= border;
	m_GridMax = box.max;
	m_GridMax += border;
	m_GridSize = m_GridMax;
	m_GridSize -= m_GridMin;

	m_GridCellSize = world_cellsize;
	//网格规格
	m_GridRes.x = (int)ceil(m_GridSize.x / world_cellsize);
	m_GridRes.y = (int)ceil(m_GridSize.y / world_cellsize);
	m_GridRes.z = (int)ceil(m_GridSize.z / world_cellsize);

	//将网格大小调整为单元格大小的倍数
	m_GridSize.x = m_GridRes.x * cell_size / sim_scale;
	m_GridSize.y = m_GridRes.y * cell_size / sim_scale;
	m_GridSize.z = m_GridRes.z * cell_size / sim_scale;
	
	//计算偏移量
	m_GridDelta = m_GridRes;
	m_GridDelta /= m_GridSize;

	int gridTotal = (int)(m_GridRes.x * m_GridRes.y * m_GridRes.z);

	rex[0] = m_GridRes.x * 8;
	rex[1] = m_GridRes.y * 8;
	rex[2] = m_GridRes.z * 8;

	m_gridData.resize(gridTotal);
}

void GridContainer::insertParticles(PointBuffer* pointBuffer) {
	std::fill(m_gridData.begin(), m_gridData.end(), -1);
	Point* p = pointBuffer->get(0);
	for (unsigned int n = 0; n < pointBuffer->size(); n++, p++) {
		int gs = getGridCellIndex(p->pos.x, p->pos.y, p->pos.z);

		//每个网格内的点划分为一个链表(m_gridData[gs]是该网格中链表的头节点)
		if (gs >= 0 && gs < m_gridData.size()) {
			p->next = m_gridData[gs];
			m_gridData[gs] = n;
		}
		else p->next = -1;
	}
}

int GridContainer::findCell(const glm::vec3& p) {
	int gc = getGridCellIndex(p.x, p.y, p.z);
	if (gc < 0 || gc >= m_gridData.size())return -1;
	return gc;
}

void GridContainer::findCells(const glm::vec3& p, float radius, int* gridCell) {
	for (int i = 0; i < 8; i++)gridCell[i] = -1;

	//计算当前粒子点光滑核所在网格范围
	int sph_min_x = ((-radius + p.x - m_GridMin.x) * m_GridDelta.x);
	int sph_min_y = ((-radius + p.y - m_GridMin.y) * m_GridDelta.y);
	int sph_min_z = ((-radius + p.z - m_GridMin.z) * m_GridDelta.z);
	if (sph_min_x < 0) sph_min_x = 0;
	if (sph_min_y < 0) sph_min_y = 0;
	if (sph_min_z < 0) sph_min_z = 0;

	//获取8个网格
	gridCell[0] = (sph_min_z * m_GridRes.y + sph_min_y) * m_GridRes.x + sph_min_x;
	gridCell[1] = gridCell[0] + 1;
	gridCell[2] = gridCell[0] + m_GridRes.x;
	gridCell[3] = gridCell[2] + 1;

	if (sph_min_z + 1 < m_GridRes.z) {
		gridCell[4] = gridCell[0] + m_GridRes.y * m_GridRes.x;
		gridCell[5] = gridCell[4] + 1;
		gridCell[6] = gridCell[4] + m_GridRes.x;
		gridCell[7] = gridCell[6] + 1;
	}
	if (sph_min_x + 1 >= m_GridRes.x) {
		gridCell[1] = -1;
		gridCell[3] = -1;
		gridCell[5] = -1;
		gridCell[7] = -1;
	}
	if (sph_min_y >= m_GridRes.y) {
		gridCell[2] = -1;
		gridCell[4] = -1;
		gridCell[6] = -1;
		gridCell[8] = -1;
	}
}

void GridContainer::findTwoCells(const glm::vec3& p, float radius, int* gridCell) {
	for (int i = 0; i < 64; i++)gridCell[i] = -1;

	//计算当前粒子点光滑核所在网格范围
	int sph_min_x = ((-radius + p.x - m_GridMin.x) * m_GridDelta.x);
	int sph_min_y = ((-radius + p.y - m_GridMin.y) * m_GridDelta.y);
	int sph_min_z = ((-radius + p.z - m_GridMin.z) * m_GridDelta.z);
	if (sph_min_x < 0) sph_min_x = 0;
	if (sph_min_y < 0) sph_min_y = 0;
	if (sph_min_z < 0) sph_min_z = 0;

	int base = (sph_min_z * m_GridRes.y + sph_min_y) * m_GridRes.x + sph_min_x;

	for (int z = 0; z < 4; z++) {
		for (int y = 0; y < 4; y++) {
			for (int x = 0; x < 4; x++) {
				if ((sph_min_x + x >= m_GridRes.x) || (sph_min_y + y >= m_GridRes.y || (sph_min_z + z >= m_GridRes.z)))
					gridCell[16 * z + 4 * y + x] = -1;
				else
					gridCell[16 * z + 4 * y + x] = base + (z * m_GridRes.y + y) * m_GridRes.x + x;
			}
		}
	}
}
