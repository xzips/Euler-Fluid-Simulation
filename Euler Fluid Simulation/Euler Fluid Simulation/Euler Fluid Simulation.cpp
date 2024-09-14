#include <iostream>
#include <vector>
#include "FluidSim.hpp"
#include <thread>
#include "SFML/Graphics.hpp"


int main()
{
	SimParameters simParams;
	DisplayParameters displayParams;


	size_t imageCount = 0;

	FluidSim fluidSim = FluidSim(simParams, displayParams);
	
	//window
	sf::RenderWindow window(sf::VideoMode(800, 800), "Euler Fluid Simulation");

	//framerate limit to 60
	window.setFramerateLimit(600);

	//create RenderTexture
	sf::Image image;
	image.create(800, 800);

	sf::Texture texture;

	texture.loadFromImage(image);


	sf::Sprite sprite(texture);

	//fluidSim.SetObstacle(16, 16, 10);



	float inVel = 4;

	for (size_t i = 0; i < simParams.n_x; i++) {
		for (size_t j = 0; j < simParams.n_y; j++) {
			float s = 1.0;	// fluid
			if (i == 0 || j == 0 || j == simParams.n_y - 1)
				s = 0.0;	// solid
			fluidSim.sField[i * simParams.n_y + j] = s;

			if (i == 1) {
				fluidSim.uField[i * simParams.n_y + j] = inVel;
			}
		}
	}



	fluidSim.SetObstacle(0.25, 0.5, 0.08);








	
	//create text object for drawing cell values
	/*
	sf::Font font;
	if (!font.loadFromFile("C:\\Users\\aspen\\Desktop\\SFMLFonts\\Helvetica.ttf"))
	{
		std::cout << "Error loading font" << std::endl;
	}
	*/
	

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();

			//if esc key is pressed, close the window
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();
		}

		//clear the window
		window.clear();








		
		float min_M = INFINITY;
		float max_M = -INFINITY;

		//find min and max
		for (size_t i = 1; i < simParams.n_x; i++)
		{
			for (size_t j = 1; j < simParams.n_y; j++)
			{
				if (fluidSim.mField[i * simParams.n_y + j] < min_M)
					min_M = fluidSim.mField[i * simParams.n_y + j];
				if (fluidSim.mField[i * simParams.n_y + j] > max_M)
					max_M = fluidSim.mField[i * simParams.n_y + j];
			}
		}


		float clipM = 0.1;

		//clip p to max 1000
		if (max_M < clipM)
			max_M = clipM;

		//for (size_t j = minJ; j < maxJ; j++)
		//	fluidSim.mField[2 * simParams.n_y + j] = 0.0;


		size_t smokeStreamThickness = 6;

		size_t j_start = simParams.n_y/2 - smokeStreamThickness;
		size_t j_end = simParams.n_y/2 + smokeStreamThickness;

		//add a rectangle of dye on the left middle
		for (size_t i = 0; i < 2; i++)
			for (size_t j = j_start; j < j_end; j++)
				fluidSim.mField[i * simParams.n_y + j] = 0.0;




		sf::Color color;

		//for each pixel in the texture, set the color based on the mapped pressure value
		for (size_t i = 0; i < window.getSize().x; i++)
		{
			for (size_t j = 0; j < window.getSize().y; j++)
			{

				size_t mapped_i = i * simParams.n_x / window.getSize().x;
				size_t mapped_j = j * simParams.n_y / window.getSize().y;
				

				

				float M_norm = 1- (fluidSim.mField[mapped_i * simParams.n_y + mapped_j] - min_M) / (max_M - min_M);


				if (M_norm < 0)
					M_norm = 0;

				if (M_norm > 1)
					M_norm = 1;
				
				color.r = (255 * M_norm);
				color.g = (255 * M_norm);
				color.b = (255 * M_norm);
			



				if (fluidSim.sField[mapped_i * simParams.n_y + mapped_j] == 0)
					color = sf::Color(28, 91, 152, 255);


				image.setPixel(i, j, color);
			}
			
		}

		//update the texture
		texture.loadFromImage(image);
		


	
		

	

		//simulate the fluid
		fluidSim.Simulate(0.005f);

		window.draw(sprite);

		//save the currnet image to up 2 folders and into "images"
		//image.saveToFile("C:\\Users\\aspen\\Desktop\\Euler-Fluid-Simulation\\Images\\image" + std::to_string(imageCount) + ".png");

		imageCount += 1;


		//draw the window
		window.display();
		



	}



}
