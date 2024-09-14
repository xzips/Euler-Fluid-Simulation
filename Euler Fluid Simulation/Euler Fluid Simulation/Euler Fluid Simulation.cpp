#include <iostream>
#include <vector>
#include "FluidSim.hpp"
#include <thread>
#include "SFML/Graphics.hpp"


int main()
{
	SimParameters simParams;
	DisplayParameters displayParams;


	FluidSim fluidSim = FluidSim(simParams, displayParams);
	
	//window
	sf::RenderWindow window(sf::VideoMode(1200, 800), "Euler Fluid Simulation");

	//framerate limit to 60
	window.setFramerateLimit(600);

	//create RenderTexture
	sf::Image image;
	image.create(1200, 800);

	sf::Texture texture;

	texture.loadFromImage(image);


	sf::Sprite sprite(texture);

	//fluidSim.SetObstacle(16, 16, 10);



	float inVel = 3;

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



	fluidSim.SetObstacle(0.35, 0.5, 0.10);




	float pipeH = 0.05 * simParams.n_y;
	size_t minJ = std::floor(0.5 * simParams.n_y - 0.5 * pipeH);
	size_t maxJ = std::floor(0.5 * simParams.n_y + 0.5 * pipeH);




	
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








		
		float min_P = INFINITY;
		float max_P = -INFINITY;

		//find min and max
		for (size_t i = 10; i < simParams.n_x; i++)
		{
			for (size_t j = 1; j < simParams.n_y; j++)
			{
				if (fluidSim.mField[i * simParams.n_y + j] < min_P)
					min_P = fluidSim.mField[i * simParams.n_y + j];
				if (fluidSim.mField[i * simParams.n_y + j] > max_P)
					max_P = fluidSim.mField[i * simParams.n_y + j];
			}
		}


		for (size_t j = minJ; j < maxJ; j++)
			fluidSim.mField[2 * simParams.n_y + j] = 0.0;





		sf::Color color;

		//for each pixel in the texture, set the color based on the mapped pressure value
		for (size_t i = 0; i < window.getSize().x; i++)
		{
			for (size_t j = 0; j < window.getSize().y; j++)
			{

				size_t mapped_i = i * simParams.n_x / window.getSize().x;
				size_t mapped_j = j * simParams.n_y / window.getSize().y;
				
				float P_norm = (fluidSim.mField[mapped_i * simParams.n_y + mapped_j] - min_P) / (max_P - min_P);

				color.r = (255 * P_norm);
				color.g = (255 * P_norm);
				color.b = (255 * P_norm);




				if (fluidSim.sField[mapped_i * simParams.n_y + mapped_j] == 0)
					color = sf::Color(0, 0, 255, 255);


				image.setPixel(i, j, color);
			}
			
		}

		//update the texture
		texture.loadFromImage(image);
		


	
		

	

		//simulate the fluid
		fluidSim.Simulate(0.01f);

		window.draw(sprite);


		//draw the window
		window.display();
		



	}



}
