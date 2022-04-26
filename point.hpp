#pragma once
#include<bits/stdc++.h>
#include<GLFW/glfw3.h>
#include<GL/glut.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>

#define ELEM_MAX 1e6
//粒子缓存类
struct Point {
	glm::vec3 pos;//点位置
	float density;//密度
	float pressure;//压力
	glm::vec3 accel;//加速度
	glm::vec3 velocity;//速度
	glm::vec3 velocity_eval;//最后速度
	glm::vec3 sum_grad_w;
	float sum_grad_w_dot;
	glm::vec3 forces;
	float correction_pressure;
	glm::vec3 correction_pressure_froce;
	glm::vec3 predicted_velocity;
	glm::vec3 predicted_position;
	float predicted_density;
	float density_error;
	float kernel_self;
	int next;//指向下一个点的索引
};

class PointBuffer {
public:
	PointBuffer();
	void reset(unsigned int capacity);//重置粒子点缓冲
	unsigned int size()const { return m_fluidCounts; }//返回粒子个数

	//获取索引为index的点
	Point* get(unsigned int index) { return m_FluidBuf + index; }
	const Point* get(unsigned int index) const { return m_FluidBuf + index; }

	Point* addPointReuse();//缓冲中加入新的点并返回
	virtual ~PointBuffer();

private:
	Point* m_FluidBuf;//粒子点缓存
	unsigned int m_fluidCounts;//粒子个数
	unsigned int m_BufCapcity;//粒子容量
};

PointBuffer::PointBuffer() :m_FluidBuf(0), m_fluidCounts(0), m_BufCapcity(0) {}
PointBuffer::~PointBuffer() {
	delete[] m_FluidBuf;
	m_FluidBuf = nullptr;
}

void PointBuffer::reset(unsigned int capacity) {
	m_BufCapcity = capacity;
	if (m_FluidBuf != 0) {//当点缓存不为空，则清空点缓存
		delete[] m_FluidBuf;
		m_FluidBuf = nullptr;
	}
	if (m_BufCapcity > 0)//给点缓存分配空间
		m_FluidBuf = new Point[m_BufCapcity]();
	m_fluidCounts = 0;//设置点数量为0
}

Point* PointBuffer::addPointReuse() {
	if (m_fluidCounts >= m_BufCapcity) {
		if (m_BufCapcity * 2 > ELEM_MAX) {
			//当超过上限时，返回最后值
			return m_FluidBuf + m_fluidCounts - 1;
		}
		m_BufCapcity *= 2;
		Point* new_data = new Point[m_BufCapcity]();
		memcpy(new_data, m_FluidBuf, m_fluidCounts * sizeof(Point));
		delete[] m_FluidBuf;
		m_FluidBuf = new_data;
	}
	//新的点
	Point* point = m_FluidBuf + m_fluidCounts++;
	point->pos = glm::vec3(0);
	point->next = 0;
	point->velocity = glm::vec3(0);
	point->velocity_eval = glm::vec3(0);
	point->pressure = 0;
	point->density = 0;
	point->accel = glm::vec3(0);

	return point;
}
