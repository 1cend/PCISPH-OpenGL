#pragma once
#define PI 3.141592f
#define minLoops 3
#define maxLoops 50
#define sizeFactor 1.0
#include"grid.hpp"
#include"neighbortbl.hpp"
#include"point.hpp"
using namespace std;
//����ϵͳ��
class FluidSystem {
public:
	FluidSystem();
	void init(unsigned int maxPointCounts,
		const glm::vec3& wallBox_min, const glm::vec3& wallBox_max,
		const glm::vec3& initFluidBox_min, const glm::vec3& initFluidBox_max,
		const glm::vec3& gravity) {
		_init(maxPointCounts, fBox3(wallBox_min, wallBox_max),
			fBox3(initFluidBox_min, initFluidBox_max), gravity);
	}
	//��ȡ��ĳߴ磨�ֽڣ�
	unsigned int getPointStride()const { return sizeof(Point); }
	//��ȡ�������
	unsigned int getPointCounts()const { return m_pointBuffer.size(); }
	//��ȡ����㻺��
	const Point* getPointBuf()const { return m_pointBuffer.get(0); }
	//�߼�֡
	void tick();
	void draw();
	~FluidSystem();
private:
	//��ʼ��ϵͳ
	void _init(unsigned int maxPointCounts, const fBox3& wallBox,
		const fBox3& initFluidBox, const glm::vec3& gravity);
	//������ʼҺ���
	void _addFluidVolume(const fBox3& fluidBox, float spacing);
	//�����ھӽṹ
	float _computeNeighbor(int);
	//�����ݶ�W
	void _computeGradWValues();
	//��������
	void _computeDensityErrorFactor();
	//��������
	void _computerExternalForces();
	glm::vec3 _boundaryForce(Point*);
	void _collisionHandling(Point*);
	//Ԥ�����
	void _predictionCorrection();
	//Ԥ���������ӵ��ٶȺ�λ��
	void _predictPositionAndVelocity(Point*);
	//����Ԥ����ܶȺ�ѹ��
	void _computePredictedDensityAndPressure(int);
	//����������ѹ��
	void _computeCorrectivePressureForce(int);
	//�����ٶȺ�λ��
	void _updatePosAndVel();
	void _boundaryHandling();
	//���ݳ�Ա
	PointBuffer m_pointBuffer;
	GridContainer m_gridContainer;
	NeighborTable m_neighborTable;

	//��λ�û�������(x,y,z)
	std::vector<float>posData;

	//SPH�⻬��
	float m_kernelPoly6;
	float m_kernelSpiky;
	float m_kernelViscosity;

	//��������
	float m_pointDistance;//�뾶
	float m_unitScale;//�ߴ絥λ
	float m_viscosity;//ճ��
	float m_restDensity;//��̬�ܶ�
	float m_pointMass;//����
	float m_smoothRadius;//�⻬�˰뾶
	float m_gasConstantK;//���峣��k
	float m_deltaTime;
	float m_densityErrorFactor;
	glm::vec3 m_gravityDir;//����ʸ��

	int m_rexSize[3];//����ߴ�

	fBox3 m_sphWallBox;

	bool _density_error_too_large;
	float _max_predicted_density;//�������ӵ����Ԥ���ܶ�
};

FluidSystem::FluidSystem() {
	m_unitScale = 0.004f;            // �ߴ絥λ
	m_viscosity = 1.0f;                // ճ��
	m_restDensity = 1000.f;            // �ܶ�
	m_pointMass = 0.0006f;            // ��������
	m_gasConstantK = 1.0f;                // �������巽�̳���
	m_smoothRadius = 0.01f;            // �⻬�˰뾶
	m_pointDistance = 0.01f;
	m_deltaTime = 0.003f;
	m_rexSize[0] = 0;
	m_rexSize[1] = 0;
	m_rexSize[2] = 0;

	//Poly6 Kernel
	m_kernelPoly6 = 315.0f / (64.0f * PI * pow(m_smoothRadius, 9));
	//Spiky Kernel
	m_kernelSpiky = -45.0f / (PI * pow(m_smoothRadius, 6));
	//Viscosity Kernel
	m_kernelViscosity = 45.0f / (PI * pow(m_smoothRadius, 6));
}

FluidSystem::~FluidSystem()
{
}
//���������еĵ�
void FluidSystem::_addFluidVolume(const fBox3& fluidBox, float spacing) {
	for (float z = fluidBox.max.z; z >= fluidBox.min.z; z -= spacing)
	{
		for (float y = fluidBox.min.y; y <= fluidBox.max.y; y += spacing)
		{
			for (float x = fluidBox.min.x; x <= fluidBox.max.x; x += spacing)
			{
				Point* p = m_pointBuffer.addPointReuse();
				p->pos = glm::vec3(x, y, z);
			}
		}
	}
}
void FluidSystem::_init(unsigned int maxPointCounts, const fBox3& wallBox, const fBox3& initFluidBox, const glm::vec3& gravity) {
	m_pointBuffer.reset(maxPointCounts);
	m_sphWallBox = wallBox;
	m_gravityDir = gravity;
	m_pointDistance = pow(m_pointMass / m_restDensity, 1.0 / 3.0);//�������Ӽ��
	//�趨�����
	_addFluidVolume(initFluidBox, m_pointDistance / m_unitScale);

	m_gridContainer.init(wallBox, m_unitScale, m_smoothRadius * 2.0, 1.0, m_rexSize);//��������ߴ�(2r)
	posData = std::vector<float>(3 * m_pointBuffer.size(), 0);
	cout << m_pointBuffer.size() << endl;
	m_gridContainer.insertParticles(&m_pointBuffer);//ÿ֡ˢ������λ��
	m_neighborTable.reset(m_pointBuffer.size());//�����ڽӱ�

	_computeGradWValues();
	_computeDensityErrorFactor();
}
float FluidSystem::_computeNeighbor(int i) {
	float h2 = m_smoothRadius * m_smoothRadius;//h^2
	Point* pi = m_pointBuffer.get(i);
	float sum = 0.0;
	m_neighborTable.point_prepare(i);
	int gridCell[8];
	m_gridContainer.findCells(pi->pos, m_smoothRadius / m_unitScale, gridCell);
	for (int cell = 0; cell < 8; cell++) {
		if (gridCell[cell] == -1) continue;
		int pndx = m_gridContainer.getGridData(gridCell[cell]);
		while (pndx != -1) {
			Point* pj = m_pointBuffer.get(pndx);
			if (pj == pi)sum += pow(h2, 3.0);
			else {
				glm::vec3 pi_pj = (pi->pos - pj->pos) * m_unitScale;
				float r2 = pi_pj.x * pi_pj.x + pi_pj.y * pi_pj.y + pi_pj.z * pi_pj.z;
				if (h2 > r2) {
					float h2_r2 = h2 - r2;
					sum += pow(h2_r2, 3.0);
					if (!m_neighborTable.point_add_neighbor(pndx, glm::sqrt(r2)))
						return sum * m_kernelPoly6;
				}
			}
			pndx = pj->next;
		}
	}
	return sum * m_kernelPoly6;
}
void FluidSystem::_computeGradWValues() {
	float h2 = m_smoothRadius * m_smoothRadius;//h^2
	const int numP = m_pointBuffer.size();
	for (int i = 0; i < numP; i++) {
		Point* pi = m_pointBuffer.get(i);
		pi->sum_grad_w = glm::vec3(0.0f);
		pi->sum_grad_w_dot = 0.0f;
	}
	for (int i = 0; i < numP; i++) {
		Point* pi = m_pointBuffer.get(i);
		//�ھ�����װ��
		pi->kernel_self = _computeNeighbor(i);
		m_neighborTable.point_commit();
		int neighborCounts = m_neighborTable.getNeighborCounts(i);
		//Ԥ���ܶȼ���
		for (int j = 0; j < neighborCounts; j++) {
			unsigned int neighborIndex;
			float r;
			m_neighborTable.getNeighborInfo(i, j, neighborIndex, r);
			Point* pj = m_pointBuffer.get(neighborIndex);
			//��Ҫ��Ԥ���λ��ȥ���¼������ͺ�ֵ
			glm::vec3 pi_pj = (pi->pos - pj->pos) * m_unitScale;
			float r2 = pi_pj.x * pi_pj.x + pi_pj.y * pi_pj.y + pi_pj.z * pi_pj.z;
			if (h2 > r2) {
				float h2_r2 = h2 - r2;
				r = pow(r2, 0.5f);
				float h = m_smoothRadius;

				glm::vec3 gradVec = (pi->pos - pj->pos) * m_kernelSpiky / r * (h - r) * (h - r);
				pi->sum_grad_w += gradVec;
				pj->sum_grad_w -= gradVec;

				pi->sum_grad_w_dot += glm::dot(gradVec, gradVec);
				pj->sum_grad_w_dot += glm::dot(-1.0f * gradVec, -1.0f * gradVec);

			}
		}
	}
}
void FluidSystem::_computeDensityErrorFactor() {
	float h2 = m_smoothRadius * m_smoothRadius;
	int maxNeighborIndex = 0;
	int maxNeighs = 0;
	const int numParticles = m_pointBuffer.size();
	for (int i = 0; i < numParticles; i++) {
		Point* pi = m_pointBuffer.get(i);
		int neighborCounts = m_neighborTable.getNeighborCounts(i);
		int numNeighbors = 0;
		//Ԥ���ܶȼ���
		for (int j = 0; j < neighborCounts; j++) {
			unsigned int neighborIndex;
			float r;
			m_neighborTable.getNeighborInfo(i, j, neighborIndex, r);
			Point* pj = m_pointBuffer.get(neighborIndex);
			//��Ҫ��Ԥ���λ��ȥ���¼������ͺ�ֵ
			glm::vec3 pi_pj = (pi->pos - pj->pos) * m_unitScale;
			float r2 = pi_pj.x * pi_pj.x + pi_pj.y * pi_pj.y + pi_pj.z * pi_pj.z;
			if (h2 > r2) {
				numNeighbors++;
			}
		}
		//��ȡ�ھ��������ӣ����ھӸ���
		if (numNeighbors > maxNeighs) {
			maxNeighs = numNeighbors;
			maxNeighborIndex = i;
		}
	}
	//��ȡ�ھ���������
	Point* pmax = m_pointBuffer.get(maxNeighborIndex);

	//�����µ�ѹ�������ܶ����
	float restVol = m_pointMass / m_restDensity;
	float preFactor = 2.0 * restVol * restVol * m_deltaTime * m_deltaTime;
	float gradWTerm = glm::dot(pmax->sum_grad_w * (-1.0f), pmax->sum_grad_w) - pmax->sum_grad_w_dot;
	float divisor = preFactor * gradWTerm;
	m_densityErrorFactor = -1.0 / divisor;
	const float factor = m_deltaTime / 0.0001f;
	float densityErrorFactorParameter = 0.05 * factor * factor;
	m_densityErrorFactor *= 1.0 * densityErrorFactorParameter;

}
void FluidSystem::_computerExternalForces() {
	float h2 = m_smoothRadius * m_smoothRadius;//h^2
	for (int i = 0; i < m_pointBuffer.size(); i++) {
		Point* pi = m_pointBuffer.get(i);
		pi->forces = glm::vec3(0.0);

		//�ھ�����װ��
		pi->kernel_self = _computeNeighbor(i);
		m_neighborTable.point_commit();

		//��������
		int neighborCounts = m_neighborTable.getNeighborCounts(i);
		const float restVolume = m_pointMass / m_restDensity;
		for (unsigned int j = 0; j < neighborCounts; j++) {
			unsigned int neighborIndex;
			float r;
			m_neighborTable.getNeighborInfo(i, j, neighborIndex, r);
			Point* pj = m_pointBuffer.get(neighborIndex);
			float h_r = m_smoothRadius - r;

			//F_Viscosity
			float vterm = m_kernelViscosity * m_viscosity * h_r * restVolume * restVolume;
			pi->forces -= (pi->velocity - pj->velocity) * vterm;
		}
		//F_gravity
		pi->forces += m_gravityDir * m_pointMass;

		//F_boundary
		pi->forces += _boundaryForce(pi) * m_pointMass;
		//��ʼ����������
		pi->correction_pressure = 0.0f;
		pi->correction_pressure_froce = glm::vec3(0.0);
	}
}
glm::vec3 FluidSystem::_boundaryForce(Point* p) {
	const float forceDistance = m_smoothRadius;
	const float invForceDist = 1.0 / forceDistance;
	const float forceStrength = sizeFactor * m_smoothRadius* m_smoothRadius*16;
	float distToWall, factor;
	glm::vec3 force(0.0);

	// ground-ceiling
	if (p->pos.y < m_sphWallBox.min.y + forceDistance) {
		distToWall = m_sphWallBox.min.y + forceDistance - p->pos.y;
		factor = distToWall * invForceDist* 2.0;
		glm::vec3 norm(0, 1.0, 0);
		force += norm * factor * forceStrength;
	}
	else if (p->pos.y > m_sphWallBox.max.y - forceDistance) {
		distToWall = p->pos.y - (m_sphWallBox.max.y - forceDistance);
		factor = distToWall * invForceDist;
		glm::vec3 norm(0, -1.0, 0);
		force +=  norm * factor * forceStrength;
	}

	// x, z boundaries
	// xy plane
	if (p->pos.x < m_sphWallBox.min.x + forceDistance) {
		distToWall = m_sphWallBox.min.x + forceDistance - p->pos.x;
		factor = distToWall * invForceDist;
		glm::vec3 norm(1.0, 0, 0);
		force += norm * factor * forceStrength;
	}
	if (p->pos.x > m_sphWallBox.max.x - forceDistance) {
		distToWall = p->pos.x - (m_sphWallBox.max.x - forceDistance);
		factor = distToWall * invForceDist;
		glm::vec3 norm(-1.0, 0, 0);
		force += norm * factor * forceStrength;   
	}

	// yz plane
	if (p->pos.z < m_sphWallBox.min.z + forceDistance) {
		distToWall = m_sphWallBox.min.z + forceDistance - p->pos.z;
		factor = distToWall * invForceDist;
		glm::vec3 norm(0, 0, 1.0);
		force += norm * factor * forceStrength;
	}
	if (p->pos.z > m_sphWallBox.max.z - forceDistance) {
		distToWall = p->pos.z - (m_sphWallBox.max.z - forceDistance);
		factor = distToWall * invForceDist;
		glm::vec3 norm(0, 0, -1.0);
		force += norm * factor * forceStrength;
	}

	return force;
}
void FluidSystem::_predictionCorrection() {
	_density_error_too_large = true;
	int iteration = 0;
	while ((iteration < minLoops) || ((_density_error_too_large) && (iteration < maxLoops))) {
		for (int i = 0; i < m_pointBuffer.size(); i++) {
			Point* p = m_pointBuffer.get(i);
			_predictPositionAndVelocity(p);
		}
		//ѭ����ֹ����
		_max_predicted_density = 0.0;

		for (int i = 0; i < m_pointBuffer.size(); i++) {
			_computePredictedDensityAndPressure(i);
		}

		//ѭ����ֹ����
		float densityErrorInPercent = max(0.1f * _max_predicted_density - 100.0f, 0.0f);
		float maxDensityErrorAllowed = 1.0f;

		//����ܶ����С����ֹ����������ֹѭ��(С���ܶ������ֵ)
		if (densityErrorInPercent < maxDensityErrorAllowed)
			_density_error_too_large = false;

		for (int i = 0; i < m_pointBuffer.size(); i++) {
			_computeCorrectivePressureForce(i);
		}
		iteration++;
	}
}
void FluidSystem::_predictPositionAndVelocity(Point* p) {
	// v' = v + delta_t * a
	// a = F / m
	// v' = v + delta_t * F * V / m
	// v' = v + delta_t * F * 1/density
	//����Ԥ���ٶȺ�λ��
	p->predicted_velocity = p->velocity + (p->forces + p->correction_pressure_froce) * m_deltaTime / m_pointMass;
	p->predicted_position = p->pos + p->predicted_velocity * m_deltaTime;
	//��ײ����
	_collisionHandling(p);
	//��ʼ��Ԥ���ܶ�
	p->predicted_density = 0.0;
	//check if particle has valid position(== not nan) to detect instabilities
	if (!(p->predicted_position[0] == p->predicted_position[0]) ||
		!(p->predicted_position[1] == p->predicted_position[1]) ||
		!(p->predicted_position[2] == p->predicted_position[2]))
	{
		std::cout << "Particle has invalid predictedPosition!!" << std::endl;
		abort();
	}
}
void FluidSystem::_collisionHandling(Point* p) {
	const float damping = 0.0;

	// collision handling with domain boundary
	for (int i = 0; i < 3; i++) {
		// check minimum box boundary
		if (p->predicted_position[i] < m_sphWallBox.min[i]) {
			p->predicted_position[i] = m_sphWallBox.min[i];
			p->predicted_velocity *= damping;
		}
		// check maximum box boundary
		if(p->predicted_position[i] > m_sphWallBox.max[i]) {
			p->predicted_position[i] = m_sphWallBox.max[i];
			p->predicted_velocity *= damping;
		}
	}
}
void FluidSystem::_computePredictedDensityAndPressure(int i) {
	float h2 = m_smoothRadius * m_smoothRadius;
	Point* pi = m_pointBuffer.get(i);
	int neighborCounts = m_neighborTable.getNeighborCounts(i);
	pi->predicted_density = pi->kernel_self * m_pointMass;
	//Ԥ���ܶȼ���
	for (int j = 0; j < neighborCounts; j++) {
		unsigned int neighborIndex;
		float r;
		m_neighborTable.getNeighborInfo(i, j, neighborIndex, r);
		Point* pj = m_pointBuffer.get(neighborIndex);
		//��Ҫ��Ԥ���λ��ȥ���¼������ͺ�ֵ
		glm::vec3 pi_pj = (pi->predicted_position - pj->predicted_position) * m_unitScale;
		float r2 = pi_pj.x * pi_pj.x + pi_pj.y * pi_pj.y + pi_pj.z * pi_pj.z;
		if (h2 > r2) {
			float h2_r2 = h2 - r2;
			pi->predicted_density += m_kernelPoly6 * pow(h2_r2, 3.0) * m_pointMass;
		}
	}
	auto sss = pi->predicted_density;
	//�����ܶ����������������
	pi->density_error = max(pi->predicted_density - m_restDensity, 0.0f);

	//����ѹ��,ֻ������ѹ����densityErrorFactor��Ԥ�ȼ��㲢��������
	pi->correction_pressure += max(pi->density_error * m_densityErrorFactor, 0.0f);

	_max_predicted_density = max(_max_predicted_density, pi->predicted_density);
}
void FluidSystem::_computeCorrectivePressureForce(int i) {
	float h2 = m_smoothRadius * m_smoothRadius;
	Point* pi = m_pointBuffer.get(i);
	int neighborCounts = m_neighborTable.getNeighborCounts(i);
	pi->correction_pressure_froce = glm::vec3(0.0f);

	float densSq = m_restDensity * m_restDensity;
	float pi_pres = pi->correction_pressure;

	//����ѹ��
	for (int j = 0; j < neighborCounts; j++) {
		unsigned int neighborIndex;
		float r;
		m_neighborTable.getNeighborInfo(i, j, neighborIndex, r);
		Point* pj = m_pointBuffer.get(neighborIndex);

		glm::vec3 pi_pj = (pi->pos - pj->pos) * m_unitScale;
		float r2 = pi_pj.x * pi_pj.x + pi_pj.y * pi_pj.y + pi_pj.z * pi_pj.z;
		if (h2 > r2) {
			r = pow(r2, 0.5);
			float h = m_smoothRadius;
			float pj_pres = pj->correction_pressure;
			glm::vec3 gradientKer = pi_pj * m_kernelSpiky / r * (h - r) * (h - r);
			float grad = pi_pres / densSq + pj_pres / (m_restDensity * m_restDensity);
			pi->correction_pressure_froce -= m_pointMass * m_pointMass * grad * gradientKer;
		}
	}
}
void FluidSystem::_updatePosAndVel() {
	for (int i = 0; i < m_pointBuffer.size(); i++) {
		Point* pi = m_pointBuffer.get(i);
		pi->forces += pi->correction_pressure_froce;
		const float invMass = 1.0 / m_pointMass;

		//Symplectic Euler Integration
		        pi->velocity+=pi->forces*invMass*m_deltaTime;
		        pi->pos+=pi->velocity*m_deltaTime/m_unitScale;

		//Leapfrog integration
		//glm::vec3 vnext = pi->velocity + pi->forces * invMass * m_deltaTime; // v(t+1/2) = v(t-1/2) + a(t) dt
		//pi->velocity_eval = (pi->velocity + vnext) * 0.5f;  //v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5
		//pi->velocity = vnext;
		//pi->pos += vnext * m_deltaTime / m_unitScale;  // p(t+1) = p(t) + v(t+1/2) dt
	}
}
void FluidSystem::_boundaryHandling(){
	const float damping = -0.1;
	for (int i = 0; i < m_pointBuffer.size(); i++)
	{
		Point* pi = m_pointBuffer.get(i);

		for (int i = 0; i < 3; i++) {
			// check minimum box boundary
			if (pi->pos[i] < m_sphWallBox.min[i]) {
				pi->pos[i] = m_sphWallBox.min[i];
				pi->velocity *= damping;
			}
			// check maximum box boundary
			if (pi->pos[i] > m_sphWallBox.max[i]) {
				pi->pos[i] = m_sphWallBox.max[i];
				pi->velocity *= damping;
			}
		}
		//����λ������
		posData[3 * i] = pi->pos.x;
		posData[3 * i + 1] = pi->pos.y;
		posData[3 * i + 2] = pi->pos.z;
	}
}
void FluidSystem::tick() {
	m_gridContainer.insertParticles(&m_pointBuffer);//ÿ֡ˢ������λ��
	m_neighborTable.reset(m_pointBuffer.size());//�����ڽӱ�

	_computerExternalForces();
	_predictionCorrection();
	_updatePosAndVel();
	_boundaryHandling();
	draw();
}
void FluidSystem::draw() {
	for (int i = 0; i < m_pointBuffer.size(); i++) {
		glPushMatrix();
		glTranslatef(posData[3 * i] * m_unitScale, posData[3 * i + 1] * m_unitScale, posData[3 * i + 2] * m_unitScale);
		glutSolidSphere(m_pointDistance, 30, 30);
		glPopMatrix();
	}

	glDisable(GL_LIGHTING);

	glColor3f(1.0f, 1.0f, 1.0f);

	// Draw bottom surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glEnd();

	// Draw top surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glEnd();

	// Draw left surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glVertex3f(m_sphWallBox.min.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glEnd();

	// Draw right surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.max.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.max.z * m_unitScale);
	glVertex3f(m_sphWallBox.max.x * m_unitScale, m_sphWallBox.min.y * m_unitScale, m_sphWallBox.min.z * m_unitScale);
	glEnd();

	glColor3f(1.0f, 0.0f, 0.0f);
	glEnable(GL_LIGHTING);
}