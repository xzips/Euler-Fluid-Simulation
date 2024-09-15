#include "FluidSim.hpp"

#include <string>

#include <iostream>



#include <execution>
#include <vector>


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


/*
* 
* This is the old, unoptimized version that takes 25ms to render a 250x120 grid onto an 800x1600 canvas
* 
void FluidSim::Render(sf::RenderWindow& window)
{

	if (displayParams.showSmoke)
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

		//for each pixel in the texture, set the color based on the mapped pressure value
		for (size_t i = 0; i < window.getSize().x; i++)
		{
			for (size_t j = 0; j < window.getSize().y; j++)
			{

				size_t mapped_i = i * simParams.n_x / window.getSize().x;
				size_t mapped_j = j * simParams.n_y / window.getSize().y;




				float M_norm = 1 - (mField[mapped_i * simParams.n_y + mapped_j] - min_M) / (max_M - min_M);


				if (M_norm < 0)
					M_norm = 0;

				if (M_norm > 1)
					M_norm = 1;

				color.r = (255 * M_norm);
				color.g = (255 * M_norm);
				color.b = (255 * M_norm);




				if (sField[mapped_i * simParams.n_y + mapped_j] == 0)
				{
					//color = sf::Color(28, 91, 152, 255);

				}


				image.setPixel(i, j, color);
			}

		}

	}


	//update the texture
	texture.loadFromImage(image);

	//draw sprite
	window.draw(sprite);

	//for obstacles DrawObstaclePretty
	for (auto& obs : simParams.obstacles)
	{
		obs.DrawObstaclePretty(window);
	}

	


	
}
*/


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

	if (displayParams.showSmoke)
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

				color.r = (255 * M_norm);
				color.g = (255 * M_norm);
				color.b = (255 * M_norm);
				

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
//takes about 450ms to solve a 600x400 grid with 200 steps
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


/*
//newer version, takes 210ms to solve a 600x400 grid with 200 steps, still too slow
void FluidSim::SolveIncompressibility(size_t numIterations, float dt) {
	size_t n_x = simParams.n_x;
	size_t n_y = simParams.n_y;
	size_t n = n_y; // Assuming 'n' is the same as 'simParams.n_y'

	float cp = simParams.density * simParams.gridSpacingX / dt;

	// Precompute red and black indices
	std::vector<size_t> redIndices;
	std::vector<size_t> blackIndices;

	// Reserve space to avoid reallocations
	redIndices.reserve((n_x - 2) * (n_y - 2) / 2);
	blackIndices.reserve((n_x - 2) * (n_y - 2) / 2);

	for (size_t i = 1; i < n_x - 1; ++i) {
		for (size_t j = 1; j < n_y - 1; ++j) {
			size_t idx = i * n + j;
			if ((i + j) % 2 == 0) {
				redIndices.push_back(idx);
			}
			else {
				blackIndices.push_back(idx);
			}
		}
	}


	//std::execution::parallel_policy exec = std::execution::par;
	std::execution::sequenced_policy exec = std::execution::seq;

	// Iterate for the specified number of iterations
	for (size_t iteration = 0; iteration < numIterations; ++iteration) {

	
		
		// --- Update Red Cells ---
		std::for_each(exec, redIndices.begin(), redIndices.end(),
			[&](size_t idx) {
				size_t i = idx / n;
				size_t j = idx % n;

				// Skip boundary cells
				if (sField[idx] == 0.0f)
					return;

				// Fetch neighboring sField values
				float sx0 = sField[(i - 1) * n + j];
				float sx1 = sField[(i + 1) * n + j];
				float sy0 = sField[i * n + (j - 1)];
				float sy1 = sField[i * n + (j + 1)];

				float s = sx0 + sx1 + sy0 + sy1;
				if (s <= 1e-6f)
					return;

				// Compute divergence
				float div = uField[(i + 1) * n + j] - uField[i * n + j]
					+ vField[i * n + (j + 1)] - vField[i * n + j];

				// Compute pressure correction
				float p = (-div / s) * simParams.overRelaxation;

				// Update pressure field
				pField[idx] += cp * p;

				// Update velocity fields
				uField[i * n + j] -= sx0 * p;
				uField[(i + 1) * n + j] += sx1 * p;
				vField[i * n + j] -= sy0 * p;
				vField[i * n + (j + 1)] += sy1 * p;
			}
		);
		sf::Clock clock;
		clock.restart();
		// --- Update Black Cells ---
		std::for_each(exec, blackIndices.begin(), blackIndices.end(),
			[&](size_t idx) {
				size_t i = idx / n;
				size_t j = idx % n;

				// Skip boundary cells
				if (sField[idx] == 0.0f)
					return;

				// Fetch neighboring sField values
				float sx0 = sField[(i - 1) * n + j];
				float sx1 = sField[(i + 1) * n + j];
				float sy0 = sField[i * n + (j - 1)];
				float sy1 = sField[i * n + (j + 1)];

				float s = sx0 + sx1 + sy0 + sy1;
				if (s <= 1e-6f)
					return;

				// Compute divergence
				float div = uField[(i + 1) * n + j] - uField[i * n + j]
					+ vField[i * n + (j + 1)] - vField[i * n + j];

				// Compute pressure correction
				float p = (-div / s) * simParams.overRelaxation;

				// Update pressure field
				pField[idx] += cp * p;

				// Update velocity fields
				uField[i * n + j] -= sx0 * p;
				uField[(i + 1) * n + j] += sx1 * p;
				vField[i * n + j] -= sy0 * p;
				vField[i * n + (j + 1)] += sy1 * p;
			}
		);


		//std::cout << "Iteration " << iteration << " took " << clock.getElapsedTime().asSeconds() *1000.f<< " ms" << std::endl;
		
	}
}
*/


/*
Problem with the above parellization function is that it maxes out my CPU usage (from 1 thread to 24 threads),
but only results in a 2x speedup, which is not worth it.

Already we use the red-black method which is good, but I would like to try multigrid and cojugate gradient method
and/or multigrid

*/




/*
Logical optimizations I will try to implemenet:
- Use newUField and newVField to store updated velocities, and rather than moving memory around,
I will swap the pointers of the uField and newUField, and vField and newVField. each iteration
(keeping track of them for any final swap), to make the code more parallelizable and cache friendly.
- possible drawback includes that adjacent cells will not be cascadingly updated meaning possibly more
iterations will be needed to reach convergence, but I think this is a good tradeoff for the possibility
of speeding up the code.


*/


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
			sf::Color color = sf::Color(255 * normalized, 0, 255 * (1 - normalized));
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

//v3 optimized version, takes 330ms for 600x400 at 200 iterations instead of old 450ms
//note this is slower than the paralell execution version, but still much better than original
//and has more room to grow, next step is to optimize using AVX2


/*
void FluidSim::SolveIncompressibility(size_t numIterations, float dt)
{
	size_t n = simParams.n_y;
	float cp = simParams.density * simParams.gridSpacingX / dt;




	//pre-computation takes about 1.8ms an saves 50ms of computation time later on with our benchmark 600x400 at 200itr example
	//from 450ms to 400ms for the full solve incompressibility function

	sf::Clock clock;
	clock.restart();
	


	
	for (size_t i = 1; i < simParams.n_x - 1; i++) {
		for (size_t j = 1; j < simParams.n_y - 1; j++) {
			float sx0 = sField[(i - 1) * n + j];
			float sx1 = sField[(i + 1) * n + j];
			float sy0 = sField[i * n + j - 1];
			float sy1 = sField[i * n + j + 1];
			float s = sx0 + sx1 + sy0 + sy1;
			
			if (s <= 1e-6)
				sWeightTemp[i * n + j] = 0.0;
			else
				sWeightTemp[i * n + j] = simParams.overRelaxation / s; //precomputing division is a big optimization

			//if s field itself is 0, set to 0
			if (sField[i * n + j] == 0.0)
			{
				sWeightTemp[i * n + j] = 0.0;
			}

			sx0Field[i * n + j] = sx0;
			sx1Field[i * n + j] = sx1;
			sy0Field[i * n + j] = sy0;
			sy1Field[i * n + j] = sy1;
		}
	}
	


	//borrow
	float* divField = newUField;
	

	std::cout << "sField values Pre-computation took " << clock.getElapsedTime().asSeconds() * 1000.f << " ms" << std::endl;
	
	clock.restart();
	
	


	//one thing I learned is that this gauss-seidel method needs to be done with rolling updates,
	//i.e. precomputing the divergence is not possible

	
	
	
	for (size_t iteration = 0; iteration < numIterations; iteration++) {

		for (size_t i = 1; i < simParams.n_x - 1; i++) {
			for (size_t j = 1; j < simParams.n_y - 1; j++) {

				size_t index = i * n + j;
				
				float swtmp = sWeightTemp[index];
				
				//if sWeightTemp is 0, continue
				if (swtmp == 0.0) [[unlikely]]
				{
						continue;
				}

					
				float div =
					uField[(i + 1) * n + j]
					- uField[index]
					+ vField[index + 1]
					- vField[index];


				
				//divField[index] = std::abs(div);



				float p = -div * swtmp;



				uField[index]        -= sx0Field[index] * p;
				uField[(i + 1) * n + j]  += sx1Field[index] * p;
				vField[index]        -= sy0Field[index] * p;
				vField[index + 1]    += sy1Field[index] * p;

				
				//only update on last iteration
				if (iteration == numIterations - 1) [[unlikely]] //branch optimijzation
				{
					pField[index] += cp * p;
				}


			}
		}
		
		//VisualizeField(divField, simParams.n_x, simParams.n_y);
		
	}





	

	std::cout << "SolveIncompressibility V3 main loop took: " << clock.getElapsedTime().asSeconds() * 1000.f << " ms" << std::endl;


	//delete



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

*/


void FluidSim::SolveIncompressibility(size_t numIterations, float dt)
{
	size_t n = simParams.n_y;
	float cp = simParams.density * simParams.gridSpacingX / dt;




	//pre-computation takes about 1.8ms an saves 50ms of computation time later on with our benchmark 600x400 at 200itr example
	//from 450ms to 400ms for the full solve incompressibility function

	sf::Clock clock;
	clock.restart();




	#pragma omp parallel for collapse(2)
	for (int i = 1; i < simParams.n_x - 1; ++i) {
		for (int j = 1; j < simParams.n_y - 1; ++j) {
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
	
	std::execution::parallel_policy exec = std::execution::par;


	for (size_t iteration = 0; iteration < numIterations; iteration++) {
		
		
		
		
		#pragma omp parallel for
		for (int r = 0; r < redIndices.size(); ++r) {
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


		#pragma omp parallel for
		for (int b = 0; b < blackIndices.size(); ++b) {
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






void FluidSim::Extrapolate() {
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
	size_t n = simParams.n_y;
	float u = (uField[i * n + j - 1] + uField[i * n + j] +
		uField[(i + 1) * n + j - 1] + uField[(i + 1) * n + j]) * 0.25;
	return u;

}

float FluidSim::avgV(size_t i, size_t j) {
	size_t n = simParams.n_y;
	float v = (vField[(i - 1) * n + j] + vField[i * n + j] +
		vField[(i - 1) * n + j + 1] + vField[i * n + j + 1]) * 0.25;
	return v;
}

/*
//old version takes ~8.0ms for 400x200 grid

void FluidSim::AdvectVel(float dt) {


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

	std::cout << "AdvectVel main loop: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	//aprox 0.3ms
	memcpy(uField, newUField, simParams.n_cells * sizeof(float));
	memcpy(vField, newVField, simParams.n_cells * sizeof(float));
}
*/


/*
//new version with STL parlellization ~1.7 ms for 400x200 grid
void FluidSim::AdvectVel(float dt) {


	//approx 0.7ms
	memcpy(newUField, uField, simParams.n_cells * sizeof(float));
	memcpy(newVField, vField, simParams.n_cells * sizeof(float));


	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5 * h;

	sf::Clock clock;

	// Assuming simParams.n_x and simParams.n_y are the sizes of the grid
	size_t totalElements = (simParams.n_x - 1) * (simParams.n_y - 1);

	// Create a range of indices to iterate over
	std::vector<size_t> indices(totalElements);
	std::iota(indices.begin(), indices.end(), 1);  // Generate indices from 1 to totalElements

	std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t idx) {
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
		});

	std::cout << "AdvectVel main loop: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	//aprox 0.3ms
	memcpy(uField, newUField, simParams.n_cells * sizeof(float));
	memcpy(vField, newVField, simParams.n_cells * sizeof(float));



}*/

//new version v3 with openMP for better performance, Xms for 400x200 grid
void FluidSim::AdvectVel(float dt) {
	// Approx 0.7ms
	memcpy(newUField, uField, simParams.n_cells * sizeof(float));
	memcpy(newVField, vField, simParams.n_cells * sizeof(float));

	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5 * h;

	sf::Clock clock;

	// Assuming simParams.n_x and simParams.n_y are the sizes of the grid
	size_t totalElements = (simParams.n_x - 1) * (simParams.n_y - 1);

#pragma omp parallel for
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

	std::cout << "AdvectVel main loop: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	// Approx 0.3ms
	memcpy(uField, newUField, simParams.n_cells * sizeof(float));
	memcpy(vField, newVField, simParams.n_cells * sizeof(float));
}




//old version using about 3.2ms to advect 400x200 grid

/*
void FluidSim::AdvectSmoke(float dt) {


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
}

*/

/*

//new version takes 0.8ms to advect 400x200 grid
void FluidSim::AdvectSmoke(float dt) {


	//memcpy
	memcpy(newMField, mField, simParams.n_cells * sizeof(float));




	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5 * h;


	
	// Compute the total number of elements
	size_t totalElements = (simParams.n_x - 2) * (simParams.n_y - 2);

	// Create a range of indices for parallelization
	std::vector<size_t> indices(totalElements);
	std::iota(indices.begin(), indices.end(), 0);  // Generate indices from 0 to totalElements - 1

	std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t idx) {
		// Map the 1D index back to 2D indices (i, j)
		size_t i = idx / (simParams.n_y - 2) + 1;
		size_t j = idx % (simParams.n_y - 2) + 1;

		// Perform the calculations if sField is not zero
		if (sField[i * n + j] != 0.0) {
			float u = (uField[i * n + j] + uField[(i + 1) * n + j]) * 0.5;
			float v = (vField[i * n + j] + vField[i * n + j + 1]) * 0.5;
			float x = i * h + h2 - dt * u;
			float y = j * h + h2 - dt * v;

			newMField[i * n + j] = SampleField(x, y, FieldType::M_FIELD);
		}
		});
	//memcpy
	memcpy(mField, newMField, simParams.n_cells * sizeof(float));
}
*/


void FluidSim::AdvectSmoke(float dt) {
	// Memcpy to copy the smoke field
	memcpy(newMField, mField, simParams.n_cells * sizeof(float));

	size_t n = simParams.n_y;
	float h = simParams.gridSpacingX;
	float h2 = 0.5 * h;

	// Compute the total number of elements
	size_t totalElements = (simParams.n_x - 2) * (simParams.n_y - 2);

#pragma omp parallel for
	for (int idx = 0; idx < totalElements; ++idx) {
		// Map the 1D index back to 2D indices (i, j)
		size_t i = idx / (simParams.n_y - 2) + 1;
		size_t j = idx % (simParams.n_y - 2) + 1;

		// Perform the calculations if sField is not zero
		if (sField[i * n + j] != 0.0) {
			float u = (uField[i * n + j] + uField[(i + 1) * n + j]) * 0.5;
			float v = (vField[i * n + j] + vField[i * n + j + 1]) * 0.5;
			float x = i * h + h2 - dt * u;
			float y = j * h + h2 - dt * v;

			newMField[i * n + j] = SampleField(x, y, FieldType::M_FIELD);
		}
	}

	// Memcpy to copy the updated smoke field back
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
	Extrapolate();
	std::cout << "Time taken to solve boundaries: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;
	
	
	clock.restart();
	AdvectVel(dt);
	std::cout << "Time taken to advect velocity: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

	
	clock.restart();
	AdvectSmoke(dt);
	std::cout << "Time taken to advect smoke: " << clock.getElapsedTime().asSeconds() * 1000.f << std::endl;

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
		
		float width = dyeSource.width * simParams.n_x;
		float height = dyeSource.height * simParams.n_y;

		for (size_t i = x - width / 2; i < x + width / 2; i++)
		{
			for (size_t j = y - height / 2; j < y + height / 2; j++)
			{
				if (i >= 1 && i < simParams.n_x - 1 && j >= 1 && j < simParams.n_y - 1)
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
	float cd = sqrt(2) * simParams.gridSpacingX;

	for (size_t i = 1; i < simParams.n_x - 2; i++) {
		for (size_t j = 1; j < simParams.n_y - 2; j++) {

			//sField[i * n + j] = 1.0;

			float dx = (i + 0.5) * simParams.gridSpacingX - x;
			float dy = (j + 0.5) * simParams.gridSpacingY - y;

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
		float cd = sqrt(2) * fluidSim->simParams.gridSpacingX;

		for (size_t i = 1; i < fluidSim->simParams.n_x - 2; i++) {
			for (size_t j = 1; j < fluidSim->simParams.n_y - 2; j++) {

				//sField[i * n + j] = 1.0;

				float dx = (i + 0.5) * fluidSim->simParams.gridSpacingY - x;
				float dy = (j + 0.5) * fluidSim->simParams.gridSpacingY - y;

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
		float cd = sqrt(2) * fluidSim->simParams.gridSpacingX;

		for (size_t i = 1; i < fluidSim->simParams.n_x - 2; i++) {
			for (size_t j = 1; j < fluidSim->simParams.n_y - 2; j++) {

				//sField[i * n + j] = 1.0;

				float dx = (i + 0.5) * fluidSim->simParams.gridSpacingY - x;
				float dy = (j + 0.5) * fluidSim->simParams.gridSpacingY - y;

				if (dx > -width / 2 && dx < width / 2 && dy > -height / 2 && dy < height / 2) {
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
				uField[i * simParams.n_y + j] = 0.0f;
				vField[i * simParams.n_y + j] = 0.0f;
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
