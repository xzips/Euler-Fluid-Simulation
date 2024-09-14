#include "FluidSim.hpp"

#include <string>

#include <iostream>

FluidSim::FluidSim(SimParameters simParams, DisplayParameters displayParams)
{

	this->simParams = simParams;
	this->displayParams = displayParams;



	uField = new float[simParams.n_cells];
	vField = new float[simParams.n_cells];

	sField = new float[simParams.n_cells];

	newUField = new float[simParams.n_cells];
	newVField = new float[simParams.n_cells];

	pField = new float[simParams.n_cells];

	mField = new float[simParams.n_cells];
	newMField = new float[simParams.n_cells];

	//memset all of them to 0.0
	std::fill(uField, uField + simParams.n_cells, 0.0f);
	std::fill(vField, vField + simParams.n_cells, 0.0f);
	std::fill(sField, sField + simParams.n_cells, 0.0f);
	std::fill(newUField, newUField + simParams.n_cells, 0.0f);
	std::fill(newVField, newVField + simParams.n_cells, 0.0f);
	std::fill(pField, pField + simParams.n_cells, 0.0f);
	std::fill(mField, mField + simParams.n_cells, 1.0f);
	std::fill(newMField, newMField + simParams.n_cells, 0.0f);




	//initialize mField to 1.0
	for (size_t i = 0; i < simParams.n_cells; i++)
	{
		mField[i] = 1.0;
	}
}



bool FluidSim::CheckFieldsExploded()
{
	std::string exploded_field;

	//if any of the fields contain NaN or Inf or -Inf etc print an error message and return true
	for (size_t i = 0; i < simParams.n_cells; i++)
	{
		if (std::isnan(uField[i]) || std::isinf(uField[i]))
		{
			exploded_field = "uField";
			break;
		}

		if (std::isnan(vField[i]) || std::isinf(vField[i]))
		{
			exploded_field = "vField";
			break;
		}

		if (std::isnan(sField[i]) || std::isinf(sField[i]))
		{
			exploded_field = "sField";
			break;
		}

		if (std::isnan(newUField[i]) || std::isinf(newUField[i]))
		{
			exploded_field = "newUField";
			break;
		}

		if (std::isnan(newVField[i]) || std::isinf(newVField[i]))
		{
			exploded_field = "newVField";
			break;
		}

		if (std::isnan(pField[i]) || std::isinf(pField[i]))
		{
			exploded_field = "pField";
			break;
		}

		if (std::isnan(mField[i]) || std::isinf(mField[i]))
		{
			exploded_field = "mField";
			break;
		}

		if (std::isnan(newMField[i]) || std::isinf(newMField[i]))
		{
			exploded_field = "newMField";
			break;
		}


	}

	if (exploded_field != "")
	{
		std::cerr << "Error: " << exploded_field << " exploded!" << std::endl;
		return true;
	}

	return false;

}


void FluidSim::ApplyGravity(float dt, float gravity) {
	size_t n = simParams.n_y;
	for (size_t i = 1; i < simParams.n_x; i++)
	{
		for (size_t j = 1; j < simParams.n_y - 1; j++)
		{
			//if the cell is not a boundary cell
			if (sField[i * n + j] != 0.0 && sField[i * n + j - 1] != 0.0)
			{
				//update velocity field
				vField[i * n + j] += gravity * dt;
			}
		}
	}
	
}


/*
If there is too much inflow to a cell, we need to redistribute
the inflow velocity equally to the neighboring cells.
*/
void FluidSim::SolveIncompressibility(size_t numIterations, float dt)
{
	size_t n = simParams.n_y;
	float cp = simParams.density * simParams.gridSpacing / dt;

	for (size_t iteration = 0; iteration < numIterations; iteration++) {

		for (size_t i = 1; i < simParams.n_x - 1; i++) {
			for (size_t j = 1; j < simParams.n_y - 1; j++) {

				//if the cell is a boundary cell
				if (sField[i * n + j] == 0.0)
					continue;

				//float s = sField[i * n + j];
				float sx0 = sField[(i - 1) * n + j];
				float sx1 = sField[(i + 1) * n + j];
				float sy0 = sField[i * n + j - 1];
				float sy1 = sField[i * n + j + 1];
				float s = sx0 + sx1 + sy0 + sy1;
				if (s <= 1e-6)
					continue;
				
				float div =
					uField[(i + 1) * n + j]
					- uField[i * n + j]
					+ vField[i * n + j + 1]
					- vField[i * n + j];

		
				
				float p = -div / s;
				p *= simParams.overRelaxation;

		

				pField[i * n + j] += cp * p;

				uField[i * n + j] -= sx0 * p;
				uField[(i + 1) * n + j] += sx1 * p;
				vField[i * n + j] -= sy0 * p;
				vField[i * n + j + 1] += sy1 * p;
			}
		}
	}
	



}


void FluidSim::Extrapolate() {
	size_t n = simParams.n_x;
	for (size_t i = 0; i < simParams.n_x; i++) {
		uField[i * n + 0] = uField[i * n + 1];
		uField[i * n + simParams.n_y - 1] = uField[i * n + simParams.n_y - 2];
	}
	for (size_t j = 0; j < simParams.n_y; j++) {
		vField[0 * n + j] = vField [1 * n + j];
		vField[(simParams.n_x - 1) * n + j] = vField[(simParams.n_x - 2) * n + j];
	}
	
}


float FluidSim::SampleField(float x, float y, FieldType field) {
	size_t n = simParams.n_y;
	float h = simParams.gridSpacing;
	float h1 = 1.0 / h;
	float h2 = 0.5 * h;
	
	x = std::max(std::min(x, (float)simParams.n_x * h), h);
	y = std::max(std::min(y, (float)simParams.n_y * h), h);

	float dx = 0.0;
	float dy = 0.0;

	float *f = nullptr;

	switch (field) {
	case FieldType::U_FIELD: f = uField; dy = h2; break;
	case FieldType::V_FIELD: f = vField; dx = h2; break;
	case FieldType::S_FIELD: f = sField; dx = h2; dy = h2; break;
	case FieldType::P_FIELD: f = pField; dx = h2; dy = h2; break;
	case FieldType::M_FIELD: f = mField; dx = h2; dy = h2; break;

	}

	if (f == nullptr) {
		return 0.0;
	}

	size_t x0 = std::min(std::floor((x - dx) * h1), (float)simParams.n_x - 1);
	float tx = ((x - dx) - x0 * h) * h1;
	size_t x1 = std::min(x0 + 1, simParams.n_x - 1);


	size_t y0 = std::min(std::floor((y - dy) * h1), (float)simParams.n_y - 1);
	float ty = ((y - dy) - y0 * h) * h1;
	size_t y1 = std::min(y0 + 1, simParams.n_y - 1);

	float sx = 1.0 - tx;
	float sy = 1.0 - ty;

	float val = sx * sy * f[x0 * n + y0] +
		tx * sy * f[x1 * n + y0] +
		tx * ty * f[x1 * n + y1] +
		sx * ty * f[x0 * n + y1];



	return val;
}


float FluidSim::avgU(size_t i, size_t j) {
	size_t n = simParams.n_x;
	float u = (uField[i * n + j - 1] + uField[i * n + j] +
		uField[(i + 1) * n + j - 1] + uField[(i + 1) * n + j]) * 0.25;
	return u;

}

float FluidSim::avgV(size_t i, size_t j) {
	size_t n = simParams.n_x;
	float v = (vField[(i - 1) * n + j] + vField[i * n + j] +
		vField[(i - 1) * n + j + 1] + vField[i * n + j + 1]) * 0.25;
	return v;
}


void FluidSim::AdvectVel(float dt) {
	//memcpy
	memcpy(newUField, uField, simParams.n_cells * sizeof(float));
	memcpy(newVField, vField, simParams.n_cells * sizeof(float));




	size_t n = simParams.n_y;
	float h = simParams.gridSpacing;
	float h2 = 0.5 * h;

	for (size_t i = 1; i < simParams.n_x; i++) {
		for (size_t j = 1; j < simParams.n_y; j++) {
			

			// u component
			if (sField[i * n + j] != 0.0 && sField[(i - 1) * n + j] != 0.0 && j <simParams.n_y - 1) {
				float x = i * h;
				float y = j * h + h2;
				float u = uField[i * n + j];
				float v = avgV(i, j);
				x = x - dt * u;
				y = y - dt * v;
				u = SampleField(x, y, FieldType::U_FIELD);
				newUField[i * n + j] = u;
			}
			// v component
			if (sField[i * n + j] != 0.0 && sField[i * n + j - 1] != 0.0 && i < simParams.n_x - 1) {
				float x = i * h + h2;
				float y = j * h;
				float u = avgU(i, j);
				float v = vField[i * n + j];
				x = x - dt * u;
				y = y - dt * v;
				v = SampleField(x, y, FieldType::V_FIELD);
			}				newVField[i * n + j] = v;

		}
	}



	//memcpy
	memcpy(uField, newUField, simParams.n_cells * sizeof(float));
	memcpy(vField, newVField, simParams.n_cells * sizeof(float));

	
}




void FluidSim::AdvectSmoke(float dt) {


	//memcpy
	memcpy(newMField, mField, simParams.n_cells * sizeof(float));




	size_t n = simParams.n_y;
	float h = simParams.gridSpacing;
	float h2 = 0.5 * h;

	for (size_t i = 1; i < simParams.n_x - 1; i++) {
		for (size_t j = 1; j < simParams.n_y - 1; j++) {

			if (sField[i * n + j] != 0.0) {
				float u = (uField[i * n + j] + uField[(i + 1) * n + j]) * 0.5;
				float v = (vField[i * n + j] + vField[i * n + j + 1]) * 0.5;
				float x = i * h + h2 - dt * u;
				float y = j * h + h2 - dt * v;

				newMField[i * n + j] = SampleField(x, y,FieldType::M_FIELD);
			}
		}
	}

	//memcpy
	memcpy(mField, newMField, simParams.n_cells * sizeof(float));
	

}


void FluidSim::Simulate(float dt) {

	//ApplyGravity(dt, simParams.gravity);
	
	
	for (size_t i = 0; i < simParams.n_cells; i++)
	{
		pField[i] = 0.0;
	}
	
	SolveIncompressibility(simParams.numIterations, dt);

	Extrapolate();
	AdvectVel(dt);
	AdvectSmoke(dt);
}




void FluidSim::SetObstacle(float x, float y, float r) {

	float vx = 0.0;
	float vy = 0.0;

	

	
	size_t n = simParams.n_x;
	float cd = sqrt(2) * simParams.gridSpacing;

	for (size_t i = 1; i < simParams.n_x - 2; i++) {
		for (size_t j = 1; j < simParams.n_y - 2; j++) {

			sField[i * n + j] = 1.0;

			float dx = (i + 0.5) * simParams.gridSpacing - x;
			float dy = (j + 0.5) * simParams.gridSpacing - y;

			if (dx * dx + dy * dy < r * r) {
				sField[i * n + j] = 0.0;
		
				mField[i * n + j] = 1.0;
				uField [i * n + j] = vx;
				uField[(i + 1) * n + j] = vx;
				vField[i * n + j] = vy;
				vField[i * n + j + 1] = vy;

				//std::cout << "i: " << i << " j: " << j << std::endl;
				
			}
		}
	}
	
}