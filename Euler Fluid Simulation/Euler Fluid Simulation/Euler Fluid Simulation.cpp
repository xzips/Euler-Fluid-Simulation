#include <iostream>
#include <vector>
#include "FluidSim.hpp"
#include <thread>
#include "SFML/Graphics.hpp"


int main()
{
	SimParameters simParams;
	DisplayParameters displayParams;

	displayParams.windowWidth = 800;
	displayParams.windowHeight = 800;
	displayParams.maxFps = 600;
	
	simParams.SetGridSize(500, 500);
	simParams.overRelaxation = 1.9;
	
	sf::RenderWindow window(sf::VideoMode(displayParams.windowWidth, displayParams.windowHeight), "Euler Fluid Simulation");
	window.setFramerateLimit(displayParams.maxFps);
	

	FluidSim fluidSim = FluidSim(simParams, displayParams);


	// Add dye source and obstacle
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
	
	fluidSim.SetObstacle(0.60, 0.5, 0.08);


	/*
	fluidSim.AddDyeSource(0.05, 0.62, 0.02, 0.01, 0);
	fluidSim.AddDyeSource(0.05, 0.58, 0.02, 0.01, 0);
	
	fluidSim.AddDyeSource(0.05, 0.54, 0.02, 0.01, 0);
	fluidSim.AddDyeSource(0.05, 0.50, 0.02, 0.01, 0);
	fluidSim.AddDyeSource(0.05, 0.46, 0.02, 0.01, 0);

	fluidSim.AddDyeSource(0.05, 0.42, 0.02, 0.01, 0);
	fluidSim.AddDyeSource(0.05, 0.38, 0.02, 0.01, 0);
	*/

	//add dye sources every 0.04
	for (size_t i = 10; i < simParams.n_x - 10; i++) {
		if (i % 10 == 0) {
			fluidSim.AddDyeSource(0.05, 0.05 + i * 0.005, 0.02, 0.001, 0);
		}
	}




	
	//create text object for drawing cell values
	
	sf::Font font;
	if (!font.loadFromFile("C:\\Users\\aspen\\Desktop\\SFMLFonts\\Helvetica.ttf"))
	{
		std::cout << "Error loading font" << std::endl;
	}

	sf::Text text;
	text.setFont(font);
	text.setCharacterSize(24);
	text.setFillColor(sf::Color::White);
	text.setPosition(10, 10);

	


	
	sf::Clock clock;
	float lastTime = 0;

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




		
		float currentTime = clock.getElapsedTime().asSeconds();
		clock.restart();
		float fps = 1.f / (currentTime);
		lastTime = currentTime;
		
		//set text string to GRID: 200x200, FPS: 60 for example
		text.setString("GRID: " + std::to_string(simParams.n_x) + "x" + std::to_string(simParams.n_y) + ", FPS: " + std::to_string(fps));
	
		
		
		
	
		fluidSim.Simulate(0.005f);

		fluidSim.Render(window);

		
		window.draw(text);

		//save the currnet image to up 2 folders and into "images"
		//image.saveToFile("C:\\Users\\aspen\\Desktop\\Euler-Fluid-Simulation\\Images\\image" + std::to_string(displayParams.frameCount) + ".png");

		displayParams.frameCount += 1;


		//draw the window
		window.display();
		



	}



}
