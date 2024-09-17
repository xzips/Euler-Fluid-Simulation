#include "FluidSim.hpp"

#include <string>

#include <iostream>



#include <math.h>

#include <vector>

//undefine _OPENMP
//#undef _OPENMP



#ifdef _OPENMP
	#include <omp.h>  // Include OpenMP header only if OpenMP is enabled
#endif

sf::RenderWindow* debugDrawWindow;
sf::Font debugDrawFont;
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

	sWeightTemp = new float[simParams.n_cells];
	sx0Field = new float[simParams.n_cells];
	sx1Field = new float[simParams.n_cells];
	sy0Field = new float[simParams.n_cells];
	sy1Field = new float[simParams.n_cells];


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

	//image.create(displayParams.windowWidth, displayParams.windowHeight); //for old render function
	image.create(simParams.n_x, simParams.n_y);
	texture.loadFromImage(image);

	sprite.setTexture(texture);


}

void FluidSim::Reset()
{
	//fill all arrays with their default values
	std::fill(uField, uField + simParams.n_cells, 0.0f);
	std::fill(vField, vField + simParams.n_cells, 0.0f);
	std::fill(sField, sField + simParams.n_cells, 1.0f);
	std::fill(newUField, newUField + simParams.n_cells, 0.0f);
	std::fill(newVField, newVField + simParams.n_cells, 0.0f);
	std::fill(pField, pField + simParams.n_cells, 0.0f);
	std::fill(mField, mField + simParams.n_cells, 1.0f);
	std::fill(newMField, newMField + simParams.n_cells, 0.0f);
	
}


void FluidSim::SetupWindTunnelBoundaries()
{
	for (size_t i = 0; i < simParams.n_x; i++) {
		for (size_t j = 0; j < simParams.n_y; j++) {
			float s = 1.0;	// fluid
			if (i == 0 || j == 0 || j == simParams.n_y - 1)
			{
				s = 0.0;	// solid
				sField[i * simParams.n_y + j] = s;
			}
			

			if (i == 1) {
				uField[i * simParams.n_y + j] = simParams.windTunnelSpeed;
			}
		}
	}
}




//New, optimized version that takes X 0.36ms instead of 25ms to render a 250x120 grid onto an 800x1600 canvas
void FluidSim::Render(sf::RenderWindow& window)
{

	//if n_x or n_y does not match the image texture size, that means grid was resized externally, so we need to
	//clear and re-create the image and texture
	if (image.getSize().x != simParams.n_x || image.getSize().y != simParams.n_y)
	{
		image.create(simParams.n_x, simParams.n_y, sf::Color::White);
		texture.loadFromImage(image);
		sprite.setTexture(texture);
	}

	if (displayParams.showDye)
	{
		float min_M = INFINITY;
		float max_M = -INFINITY;

		//find min and max
		for (size_t i = 1; i < simParams.n_x; i++)
		{
			for (size_t j = 1; j < simParams.n_y; j++)
			{
				if (mField[i * simParams.n_y + j] < min_M)
					min_M = mField[i * simParams.n_y + j];
				if (mField[i * simParams.n_y + j] > max_M)
					max_M = mField[i * simParams.n_y + j];
			}
		}



		sf::Color color;

		//for each sim box in the texture, set the color based on the m value
		for (size_t i = 0; i < simParams.n_x; i++)
		{
			for (size_t j = 0; j < simParams.n_y; j++)
			{
				
				
				float M_norm = 1 - (mField[i * simParams.n_y + j] - min_M) / (max_M - min_M);

				if (M_norm < 0)
					M_norm = 0;

				if (M_norm > 1)
					M_norm = 1;

				color.r = (sf::Uint8)(255 * M_norm);
				color.g = (sf::Uint8)(255 * M_norm);
				color.b = (sf::Uint8)(255 * M_norm);
				

				if (sField[i * simParams.n_y + j] == 0)
				{
					//color = sf::Color(28, 91, 152, 255);
				}


				image.setPixel(i, j, color);
			}

		}

	}


	//update the texture
	texture.loadFromImage(image);

	//rescale the sprite to fit the window
	sprite.setScale((float)window.getSize().x / (float)texture.getSize().x, (float)window.getSize().y / (float)texture.getSize().y);

	//print the window size and texture size

	//draw sprite
	window.draw(sprite);

	//for obstacles DrawObstaclePretty
	for (auto& obs : simParams.obstacles)
	{
		obs.DrawObstaclePretty(window);
	}
}







void VisualizeField(float*field, size_t n_x, size_t n_y)
{


	static sf::Image image;
	image.create(n_x, n_y, sf::Color::Black);
	

	//get min and max
	float min = std::numeric_limits<float>::max();
	float max = std::numeric_limits<float>::min();
	
	for (size_t i = 1; i < n_x -1; i++)
	{
		for (size_t j = 1; j < n_y -1; j++)
		{
			float val = field[i * n_y + j];
			if (val < min)
				min = val;
			if (val > max)
				max = val;
		}
	}


	for (size_t i = 0; i < n_x; i++)
	{
		for (size_t j = 0; j < n_y; j++)
		{
			float val = field[i * n_y + j];
			float normalized = (val - min) / (max - min);
			sf::Color color = sf::Color((sf::Uint8)(255.f * normalized), 0, (sf::Uint8)(255.f * (1 - normalized)));
			image.setPixel(i, j, color);
		}
	}

	static sf::Texture texture;
	texture.loadFromImage(image);
	static sf::Sprite sprite(texture);


	//rescale sprite to fit window
	sprite.setScale((float)debugDrawWindow->getSize().x / n_x, (float)debugDrawWindow->getSize().y / n_y);

	//create text
	sf::Text text;
	text.setFont(debugDrawFont);

	//set to the top left corner
	text.setPosition(0, 0);
	text.setCharacterSize(20);
	text.setFillColor(sf::Color::White);
	text.setString("Max Value: " + std::to_string(max));

	
	debugDrawWindow->draw(sprite);

	debugDrawWindow->draw(text);

	debugDrawWindow->display();
}

bool FluidSim::PositionIsBoundary(float x, float y)
{
	//check if a position is a boundary cell (i.e. s = 0)
	size_t n = simParams.n_y;
	size_t i = std::min(std::max((size_t)(x / simParams.gridSpacingY), (size_t)0), simParams.n_x - 1);
	size_t j = std::min(std::max((size_t)(y / simParams.gridSpacingY), (size_t)0), simParams.n_y - 1);

	return sField[i * n + j] == 0.0;




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
//takes about 450ms to solve a 600x400 grid with 200 steps, original function kept commented for readability
void FluidSim::SolveIncompressibility(size_t numIterations, float dt)
{
	size_t n = simParams.n_y;
	float cp = simParams.density * simParams.gridSpacingX / dt;



	
	for (size_t iteration = 0; iteration < numIterations; iteration++) {


		sf::Clock clock;

		clock.restart();

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
*/






#ifndef CXXDROID_COMPAT

//takes about 50ms to solve a 600x400 grid with 200 steps on a 12-core 3.8GHz Ryzen. that's approx 9x performance gain
void FluidSim::SolveIncompressibility(size_t numIterations, float dt)
{
	size_t n = simParams.n_y;
	float cp = simParams.density * simParams.gridSpacingX / dt;



	sf::Clock clock;
	clock.restart();



#ifdef _OPENMP
	#pragma omp parallel // for collapse(2) //seemed to generate a warning it was being ignored, so I removed it lol
#endif
	for (int i = 1; i < (int)simParams.n_x - 1; ++i) {
		for (int j = 1; j < (int)simParams.n_y - 1; ++j) {
			size_t idx = i * n + j;

			// Access neighboring values in the sField
			float sx0 = sField[(i - 1) * n + j];
			float sx1 = sField[(i + 1) * n + j];
			float sy0 = sField[i * n + j - 1];
			float sy1 = sField[i * n + j + 1];

			// Sum up the neighboring values
			float s = sx0 + sx1 + sy0 + sy1;

			// Precompute sWeightTemp
			if (s <= 1e-6) {
				sWeightTemp[idx] = 0.0;
			}
			else {
				sWeightTemp[idx] = simParams.overRelaxation / s; // precomputing division
			}

			// If the current sField is zero, set sWeightTemp to zero
			if (sField[idx] == 0.0) {
				sWeightTemp[idx] = 0.0;
			}

			// Update the neighboring values for the next calculations
			sx0Field[idx] = sx0;
			sx1Field[idx] = sx1;
			sy0Field[idx] = sy0;
			sy1Field[idx] = sy1;
		}
	}




	//borrow
	float* divField = newUField;


	std::cout << "sField values Pre-computation took " << clock.getElapsedTime().asSeconds() * 1000.f << " ms" << std::endl;

	clock.restart();




	//one thing I learned is that this gauss-seidel method needs to be done with rolling updates,
	//i.e. precomputing the divergence is not possible


	// Precompute red and black indices
	std::vector<size_t> redIndices;
	std::vector<size_t> blackIndices;

	
	redIndices.reserve((simParams.n_x - 2) * (simParams.n_y - 2) / 2);
	blackIndices.reserve((simParams.n_x - 2) * (simParams.n_y - 2) / 2);



	
	for (size_t i = 1; i < simParams.n_x - 1; ++i) {
		for (size_t j = 1; j < simParams.n_y - 1; ++j) {
			size_t idx = i * n + j;
			if ((i + j) % 2 == 0) {
				redIndices.push_back(idx);
			}
			else {
				blackIndices.push_back(idx);
			}
		}
	}
	


	for (size_t iteration = 0; iteration < numIterations; iteration++) {
		
		
		
#ifdef _OPENMP
		#pragma omp parallel for
#endif
		for (int r = 0; r < (int)redIndices.size(); ++r) {
			size_t idx = redIndices[r];
			size_t i = idx / n;
			size_t j = idx % n;

			float swtmp = sWeightTemp[idx];

			// If sWeightTemp is 0, continue
			if (swtmp == 0.0) {
				continue;
			}

			// Calculate divergence
			float div = uField[(i + 1) * n + j] - uField[idx] +
				vField[idx + 1] - vField[idx];

			// Compute pressure correction
			float p = -div * swtmp;

			// Update velocity fields
			uField[idx] -= sx0Field[idx] * p;
			uField[(i + 1) * n + j] += sx1Field[idx] * p;
			vField[idx] -= sy0Field[idx] * p;
			vField[idx + 1] += sy1Field[idx] * p;

			// Update pressure field only on the last iteration
			if (iteration == numIterations - 1) {
				pField[idx] += cp * p;
			}
		}

#ifdef _OPENMP
		#pragma omp parallel for
#endif
		for (int b = 0; b < (int)blackIndices.size(); ++b) {
			size_t idx = blackIndices[b];
			size_t i = idx / n;
			size_t j = idx % n;

			float swtmp = sWeightTemp[idx];

			// If sWeightTemp is 0, continue
			if (swtmp == 0.0) {
				continue;
			}

			// Calculate divergence
			float div = uField[(i + 1) * n + j] - uField[idx] +
				vField[idx + 1] - vField[idx];

			// Compute pressure correction
			float p = -div * swtmp;

			// Update velocity fields
			uField[idx] -= sx0Field[idx] * p;
			uField[(i + 1) * n + j] += sx1Field[idx] * p;
			vField[idx] -= sy0Field[idx] * p;
			vField[idx + 1] += sy1Field[idx] * p;

			// Update pressure field only on the last iteration
			if (iteration == numIterations - 1) {
				pField[idx] += cp * p;
			}
		}
		
		
	}


	std::cout << "SolveIncompressibility V4 main loop took: " << clock.getElapsedTime().asSeconds() * 1000.f << " ms" << std::endl;


	//get min and max
	float min = std::numeric_limits<float>::max();
	float max = std::numeric_limits<float>::min();

	float sum = 0.0;

	for (size_t i = 1; i < simParams.n_x - 1; i++)
	{
		for (size_t j = 1; j < simParams.n_y - 1; j++)
		{
			float val = divField[i * simParams.n_y + j];
			if (val < min)
				min = val;
			if (val > max)
				max = val;

			sum += val;
		}
	}

	simStats.maxVelocityFieldDivergence = max;
	simStats.avgVelocityFieldDivergence = sum / (simParams.n_x * simParams.n_y);


}

#endif


//CXXDroid does not have OpenMP, but we can use C++20 paralell exection for a ~2x performance improvement
#ifdef CXXDROID_COMPAT

void FluidSim::SolveIncompressibility(size_t numIterations, float dt)
{
	size_t n = simParams.n_y;
	float cp = simParams.density * simParams.gridSpacingX / dt;




	for (size_t iteration = 0; iteration < numIterations; iteration++) {


		sf::Clock clock;

		clock.restart();

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


#endif


//this function only takes 0.6ms for a 600x400 grid and turns out to be slower when paralellized, so will leave it as sequential execution
void FluidSim::SolveBoundaries() {
	size_t n = simParams.n_y;
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
	float h = simParams.gridSpacingX;
	float h1 = 1.f / h;
	float h2 = 0.5f * h;
	
	x = std::max(std::min(x, (float)simParams.n_x * h), h);
	y = std::max(std::min(y, (float)simParams.n_y * h), h);

	float dx = 0.0f;
	float dy = 0.0f;

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

	size_t x0 = (size_t)std::min(floor((x - dx) * h1), (float)simParams.n_x - 1);
	
	float tx = ((x - dx) - x0 * h) * h1;
	size_t x1 = std::min(x0 + 1, simParams.n_x - 1);


	size_t y0 = (size_t)std::min(floor((y - dy) * h1), (float)simParams.n_y - 1);
	
	
	float ty = ((y - dy) - y0 * h) * h1;
	size_t y1 = std::min(y0 + 1, simParams.n_y - 1);

	float sx = 1.0f - tx;
	float sy = 1.0f - ty;

	float val = sx * sy * f[x0 * n + y0] +
		tx * sy * f[x1 * n + y0] +
		tx * ty * f[x1 * n + y1] +
		sx * ty * f[x0 * n + y1];



	return val;
}


float FluidSim::avgU(size_t i, size_t j) {
	size_t n = simParams.n_y;
	float u = (uField[i * n + j - 1] + uField[i * n + j] +
		uField[(i + 1) * n + j - 1] + uField[(i + 1) * n + j]) * 0.25f;
	return u;

}

float FluidSim::avgV(size_t i, size_t j) {
	size_t n = simParams.n_y;
	float v = (vField[(i - 1) * n + j] + vField[i * n + j] +
		vField[(i - 1) * n + j + 1] + vField[i * n + j + 1]) * 0.25f;
	return v;
}

/*
//old version takes ~29ms for 600x400 grid
void FluidSim::AdvectVelocity(float dt) {


	//approx 0.7ms
	memcpy(newUField, uField, simParams.n_cells * sizeof(float));
	memcpy(newVField, vField, simParams.n_cells * sizeof(float));


	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5 * h;

	sf::Clock clock;

	//main loop approx 7.0ms
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
				newVField[i * n + j] = v;
			}
		}
	}

	std::cout << "AdvectVelocity main loop: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	//aprox 0.3ms
	memcpy(uField, newUField, simParams.n_cells * sizeof(float));
	memcpy(vField, newVField, simParams.n_cells * sizeof(float));
}
*/



//new version with openMP for better performance, 3.1ms for 600x400 grid, 10x peformance gain
void FluidSim::AdvectVelocity(float dt) {
	// Approx 0.7ms
	memcpy(newUField, uField, simParams.n_cells * sizeof(float));
	memcpy(newVField, vField, simParams.n_cells * sizeof(float));

	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5f * h;

	sf::Clock clock;

	// Assuming simParams.n_x and simParams.n_y are the sizes of the grid
	int totalElements = (int)(simParams.n_x - 1) * (simParams.n_y - 1);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int idx = 1; idx < totalElements; ++idx) {
		size_t i = idx / (simParams.n_y - 1) + 1;
		size_t j = idx % (simParams.n_y - 1) + 1;

		// u component
		if (sField[i * n + j] != 0.0 && sField[(i - 1) * n + j] != 0.0 && j < simParams.n_y - 1) {
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
			newVField[i * n + j] = v;
		}
	}

	std::cout << "AdvectVelocity main loop: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	// Approx 0.3ms
	memcpy(uField, newUField, simParams.n_cells * sizeof(float));
	memcpy(vField, newVField, simParams.n_cells * sizeof(float));
}




//old version using about 13.1ms to advect 600x400 grid

/*
void FluidSim::AdvectDye(float dt) {


	//memcpy
	memcpy(newMField, mField, simParams.n_cells * sizeof(float));




	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
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
}*/




//new version takes about 1.6ms to advect 600x400 grid, that's an 8x improvement
void FluidSim::AdvectDye(float dt) {
	// Memcpy to copy the dye field
	memcpy(newMField, mField, simParams.n_cells * sizeof(float));

	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5f * h;

	// Compute the total number of elements
	int totalElements = (int)(simParams.n_x - 2) * (simParams.n_y - 2);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int idx = 0; idx < totalElements; ++idx) {
		// Map the 1D index back to 2D indices (i, j)
		size_t i = idx / (simParams.n_y - 2) + 1;
		size_t j = idx % (simParams.n_y - 2) + 1;

		// Perform the calculations if sField is not zero
		if (sField[i * n + j] != 0.0) {
			float u = (uField[i * n + j] + uField[(i + 1) * n + j]) * 0.5f;
			float v = (vField[i * n + j] + vField[i * n + j + 1]) * 0.5f;
			float x = i * h + h2 - dt * u;
			float y = j * h + h2 - dt * v;

			newMField[i * n + j] = SampleField(x, y, FieldType::M_FIELD);
		}
	}

	// Memcpy to copy the updated dye field back
	memcpy(mField, newMField, simParams.n_cells * sizeof(float));
}






void FluidSim::Simulate(float dt) {
	ApplyDyeSources();

	//ApplyGravity(dt, simParams.gravity);
	
	
	for (size_t i = 0; i < simParams.n_cells; i++)
	{
		pField[i] = 0.0;
	}
	

	static sf::Clock clock;

	clock.restart();
	SolveIncompressibility(simParams.numIterations, dt);
	std::cout << "Time taken to solve incompressibility: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;
	
	clock.restart();
	SolveBoundaries();
	std::cout << "Time taken to solve boundaries: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;
	
	
	clock.restart();
	AdvectVelocity(dt);
	std::cout << "Time taken to advect velocity: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	
	clock.restart();
	AdvectDye(dt);
	std::cout << "Time taken to advect dye: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

}

void FluidSim::ApplyDyeSources()
{
	for (auto& dyeSource : simParams.dyeSources)
	{
		//for each one, assume that x and y are in world coordinates i.e. 0 to 1, and centered where rect should be
		//so need to first find the cell that the center of the rect is in, and then use two for loops to add dye to the cells in the rect
		//ensure that we check for boundary conditions
		float x = dyeSource.x * simParams.n_x;
		float y = dyeSource.y * simParams.n_y;
		
		float width = (float)(dyeSource.width * simParams.n_x);
		float height = (float)(dyeSource.height * simParams.n_y);

		for (int i = (int)(x - width / 2); i < x + width / 2; i++)
		{
			for (int j = (int)(y - height / 2); j < y + height / 2; j++)
			{
				if (i >= 1 && i < (int)simParams.n_x - 1 && j >= 1 && j < (int)simParams.n_y - 1)
				{
					mField[i * simParams.n_y + j] = dyeSource.density;
				}
			}
		}

	}
}

void FluidSim::AddDyeSource(float x, float y, float width, float height, float density = 0)
{
	RectangularDyeSource source;
	source.x = x;
	source.y = y;
	source.width = width;
	source.height = height;
	source.density = density;

	simParams.dyeSources.push_back(source);

}


void FluidSim::SetObstacle(float x, float y, float r) {

	float vx = 0.0;
	float vy = 0.0;

	

	
	size_t n = simParams.n_y;
	float cd = sqrt(2.f) * simParams.gridSpacingX;

	for (size_t i = 1; i < simParams.n_x - 2; i++) {
		for (size_t j = 1; j < simParams.n_y - 2; j++) {

			//sField[i * n + j] = 1.0;

			float dx = (i + 0.5f) * simParams.gridSpacingX - x;
			float dy = (j + 0.5f) * simParams.gridSpacingY - y;

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



void Obstacle::SetObstacleSField(FluidSim* fluidSim)
{
	if (type == ObstacleType::CIRCLE)
	{
		float vx = 0.0;
		float vy = 0.0;

		size_t n = fluidSim->simParams.n_y;
		float cd = sqrt(2.f) * fluidSim->simParams.gridSpacingX;

		for (size_t i = 1; i < fluidSim->simParams.n_x - 2; i++) {
			for (size_t j = 1; j < fluidSim->simParams.n_y - 2; j++) {

				//sField[i * n + j] = 1.0;

				float dx = (i + 0.5f) * fluidSim->simParams.gridSpacingY - x;
				float dy = (j + 0.5f) * fluidSim->simParams.gridSpacingY - y;

				if (dx * dx + dy * dy < radius * radius) {
					fluidSim->sField[i * n + j] = 0.0;

					fluidSim->mField[i * n + j] = 1.0;
					fluidSim->uField[i * n + j] = vx;
					fluidSim->uField[(i + 1) * n + j] = vx;
					fluidSim->vField[i * n + j] = vy;
					fluidSim->vField[i * n + j + 1] = vy;

					//std::cout << "i: " << i << " j: " << j << std::endl;

				}
			}
		}
	}

	if (type == ObstacleType::SQUARE)
	{
		float vx = 0.0;
		float vy = 0.0;

		size_t n = fluidSim->simParams.n_y;
		float cd = sqrt(2.f) * fluidSim->simParams.gridSpacingX;

		for (size_t i = 1; i < fluidSim->simParams.n_x - 2; i++) {
			for (size_t j = 1; j < fluidSim->simParams.n_y - 2; j++) {

				//sField[i * n + j] = 1.0;

				float dx = (i + 0.5f) * fluidSim->simParams.gridSpacingY - x;
				float dy = (j + 0.5f) * fluidSim->simParams.gridSpacingY - y;

				if (dx > -width / 2 && dx < width / 2 && dy > -height / 2 && dy < height / 2) {
					fluidSim->sField[i * n + j] = 0.0f;

					fluidSim->mField[i * n + j] = 1.0f;

					//if we are on the boundary set u and v to vx and vy
					if (dx == -width / 2 || dx == width / 2 || dy == -height / 2 || dy == height / 2)
					{
						fluidSim->uField[i * n + j] = fluidSim->simParams.mouseVelocity.x;
						fluidSim->uField[(i + 1) * n + j] = fluidSim->simParams.mouseVelocity.x;
						fluidSim->vField[i * n + j] = fluidSim->simParams.mouseVelocity.y;
						fluidSim->vField[i * n + j + 1] = fluidSim->simParams.mouseVelocity.y;
					}
					else {
						fluidSim->uField[i * n + j] = 0.0;
						fluidSim->uField[(i + 1) * n + j] = 0.0;
						fluidSim->vField[i * n + j] = 0.0;
						fluidSim->vField[i * n + j + 1] = 0.0;
					}

					//std::cout << "i: " << i << " j: " << j << std::endl;

				}
			}
		}
		
	}


	if (type == ObstacleType::IMAGE)
	{
		// Black is obstacle, white is background, use brightness 128 as threshold

		size_t n_x = fluidSim->simParams.n_x;
		size_t n_y = fluidSim->simParams.n_y;

		// Image dimensions
		size_t imgWidth = modelImage.getSize().x;
		size_t imgHeight = modelImage.getSize().y;

		// Determine scaling factors to fit the image into simulation bounds
		float simMinDim = static_cast<float>(std::min(n_x, n_y));
		float imgMaxDim = static_cast<float>(std::max(imgWidth, imgHeight));

		// Calculate the scaling factor to fit the image within the simulation bounds
		float fitScale = simMinDim / imgMaxDim;

		// Apply the regular scaling factor (local variable)
		float totalScale = fitScale * scale;

		// Adjusted image width and height after scaling
		float scaledImgWidth = imgWidth * totalScale;
		float scaledImgHeight = imgHeight * totalScale;

		// Simulation grid dimensions
		float simGridWidth = static_cast<float>(n_x);
		float simGridHeight = static_cast<float>(n_y);

		float xposInGrid = x * simGridWidth / ((float)n_x / (float)n_y);
		float yposInGrid = y * simGridHeight;

		// Compute offsets to center the image in "x" and "y" local variables
		float simOffsetX = xposInGrid - scaledImgWidth / 2;
		float simOffsetY = yposInGrid - scaledImgHeight / 2;


		for (size_t i = 1; i < n_x - 2; i++) {
			for (size_t j = 1; j < n_y - 2; j++) {

				// Compute the position relative to the image
				float simX = static_cast<float>(i) + 0.5f;
				float simY = static_cast<float>(j) + 0.5f;

				// Adjust for offset to center the image
				float imgXf = (simX - simOffsetX) / totalScale;
				float imgYf = (simY - simOffsetY) / totalScale;
				

				// Ensure we don't access out-of-bounds pixels
				if (imgXf >= 0 && imgXf < imgWidth && imgYf >= 0 && imgYf < imgHeight) {

					size_t imgX = static_cast<size_t>(imgXf);
					size_t imgY = static_cast<size_t>(imgYf);

					// Check pixel brightness threshold, alpha must be > 128
					if (modelImage.getPixel(imgX, imgY).a > 128) {
						size_t idx = i * n_y + j;
						fluidSim->sField[idx] = 0.0f;

						fluidSim->mField[idx] = 1.0f;
						fluidSim->uField[idx] = 0.0f;
						fluidSim->uField[(i + 1) * n_y + j] = 0.0f;
						fluidSim->vField[idx] = 0.0f;
						fluidSim->vField[i * n_y + j + 1] = 0.0f;
					}
				}
			}
		}
	}

		

}
	

	

void FluidSim::UpdateSField()
{
	//fill s field with 1 
	std::fill(sField, sField + simParams.n_cells, 1.0f);

	
	for (auto obstacle : simParams.obstacles)
	{
		obstacle.SetObstacleSField(this);
	}
	//anywhere that s is now 0, set u and v to 0
	for (size_t i = 0; i < simParams.n_x; i++) {
		for (size_t j = 0; j < simParams.n_y; j++) {
			if (sField[i * simParams.n_y + j] == 0.0f) {
				uField[i * simParams.n_y + j] = simParams.mouseVelocity.x;
				vField[i * simParams.n_y + j] = simParams.mouseVelocity.y;
			}
		}
	}


	SetupWindTunnelBoundaries();
	

}









void Obstacle::DrawObstaclePretty(sf::RenderWindow& window)
{


	// Get the window size
	sf::Vector2u windowSize = window.getSize();
	float windowWidth = static_cast<float>(windowSize.x);
	float windowHeight = static_cast<float>(windowSize.y);

	float domainHeight = 1.0f;
	float domainWidth = (float)windowWidth / (float)windowHeight;;


	// Compute scaling factors to map simulation units to window pixels
	float scaleX = windowWidth / domainWidth;
	float scaleY = windowHeight / domainHeight;

	if (type == ObstacleType::SQUARE)
	{
		// Convert size to window units
		float rectWidth = width * scaleX;
		float rectHeight = height * scaleY;

		// Position rectangle centered at (x, y)
		float rectX = x * scaleX - rectWidth / 2.0f;
		float rectY = y * scaleY - rectHeight / 2.0f;

		sf::RectangleShape rectangle(sf::Vector2f(rectWidth, rectHeight));
		rectangle.setPosition(rectX, rectY);

		// Set appearance
		rectangle.setFillColor(sf::Color(63, 72, 204)); // Example color
		rectangle.setOutlineColor(sf::Color::White);
		rectangle.setOutlineThickness(1.0f);

		window.draw(rectangle);
	}
	else if (type == ObstacleType::CIRCLE)
	{
		// Convert radius to window units
		float circleRadiusX = radius * scaleX;
		float circleRadiusY = radius * scaleY;
		// Use average radius for uniform scaling
		float circleRadius = (circleRadiusX + circleRadiusY) / 2.0f;

		// Position circle centered at (x, y)
		float circleX = x * scaleX - circleRadius;
		float circleY = y * scaleY - circleRadius;

		sf::CircleShape circle(circleRadius);
		circle.setPosition(circleX, circleY);

		// Set appearance
		circle.setFillColor(sf::Color(63, 72, 204)); // Example color
		circle.setOutlineColor(sf::Color::White);
		circle.setOutlineThickness(1.0f);

		window.draw(circle);
	}
	else if (type == ObstacleType::IMAGE)
	{
		modelTexture.loadFromImage(modelImage);
		sprite.setTexture(modelTexture);

		// Get image dimensions
		float imgWidth = static_cast<float>(modelImage.getSize().x);
		float imgHeight = static_cast<float>(modelImage.getSize().y);

		// Image aspect ratio
		float imgAspectRatio = imgWidth / imgHeight;

		// Desired width and height in simulation units
		// Assume 'scale' represents the desired width in simulation units
		float desiredWidthSim = scale;
		float desiredHeightSim = desiredWidthSim / imgAspectRatio;

		// Desired width and height in window units
		float desiredWidthWin = desiredWidthSim * scaleX;
		float desiredHeightWin = desiredHeightSim * scaleY;

		// Scaling factors for the sprite
		float spriteScaleX = desiredWidthWin / imgWidth;
		float spriteScaleY = desiredHeightWin / imgHeight;

		sprite.setScale(spriteScaleX, spriteScaleY);

		// Set the origin to the center of the sprite
		sprite.setOrigin(imgWidth / 2.0f, imgHeight / 2.0f);

		// Position the sprite
		float spriteX = x * scaleX;
		float spriteY = y * scaleY;
		sprite.setPosition(spriteX, spriteY);

		window.draw(sprite);
	}
	
	
}
