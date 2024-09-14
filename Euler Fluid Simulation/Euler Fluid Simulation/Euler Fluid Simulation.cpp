#include <iostream>
#include <vector>
#include "FluidSim.hpp"
#include <thread>
#include "SFML/Graphics.hpp"
#include <windows.h>
#include <string>
#include "ModelLoading.hpp"



int main()
{
	SimParameters simParams;
	DisplayParameters displayParams;

	displayParams.windowWidth = 800;
	displayParams.windowHeight = 800;
	displayParams.maxFps = 600;

	simParams.windTunnelSpeed = 4.f;
	
	simParams.SetGridSize(250, 250);
	simParams.overRelaxation = 1.9;
	
	sf::RenderWindow window(sf::VideoMode(displayParams.windowWidth, displayParams.windowHeight), "Euler Fluid Simulation");
	window.setFramerateLimit(displayParams.maxFps);
	

	FluidSim fluidSim = FluidSim(simParams, displayParams);

	

	
	
	Obstacle circle = Obstacle(0.25, 0.5, 0.08);
	//Obstacle rect = Obstacle(0.25, 0.48, 0.1, 0.1);

	fluidSim.simParams.obstacles.push_back(circle);
	fluidSim.UpdateSField();

	


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

	sf::Vector2f loadModelButtonPosition(displayParams.windowWidth - 100, 50);

	sf::RectangleShape loadModelButton(sf::Vector2f(150, 40));
	loadModelButton.setFillColor(sf::Color(70, 70, 70));
	//center origin
	loadModelButton.setOrigin(loadModelButton.getLocalBounds().width / 2, loadModelButton.getLocalBounds().height / 2);

	loadModelButton.setPosition(loadModelButtonPosition);


	sf::Text loadModelText;
	loadModelText.setFont(font);
	loadModelText.setCharacterSize(24);
	loadModelText.setFillColor(sf::Color::White);
	loadModelText.setString("Load Model");
	//center the text
	loadModelText.setOrigin(loadModelText.getLocalBounds().width / 2, loadModelText.getLocalBounds().height / 2 + 5);
	loadModelText.setPosition(loadModelButtonPosition.x, loadModelButtonPosition.y);


	
	sf::Vector2f resetButtonPosition(displayParams.windowWidth - 100, 100);
	
	sf::RectangleShape resetButton(sf::Vector2f(150, 40));
	resetButton.setFillColor(sf::Color(70, 70, 70));
	//center origin
	resetButton.setOrigin(resetButton.getLocalBounds().width / 2, resetButton.getLocalBounds().height / 2);
	
	resetButton.setPosition(resetButtonPosition);
	
	sf::Text resetText;
	resetText.setFont(font);
	resetText.setCharacterSize(24);
	resetText.setFillColor(sf::Color::White);
	resetText.setString("Reset");
	//center the text
	resetText.setOrigin(resetText.getLocalBounds().width / 2, resetText.getLocalBounds().height / 2 + 5);
	resetText.setPosition(resetButtonPosition.x, resetButtonPosition.y);
	
	

	



	sf::Clock clock;
	float lastTime = 0;

	while (window.isOpen())
	{
		
		bool leftMouseClickedThisFrame = false;
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();

			//if esc key is pressed, close the window
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();

			//if left mouse button is clicked, set leftMouseClickedThisFrame to true
			if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left)
				leftMouseClickedThisFrame = true;
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


		//draw the load model button
		window.draw(loadModelButton);
		window.draw(loadModelText);

		//draw the reset button
		window.draw(resetButton);
		window.draw(resetText);

		//if the mouse is over the button, change the color
		if (sf::Mouse::getPosition(window).x > loadModelButtonPosition.x - 75 && sf::Mouse::getPosition(window).x < loadModelButtonPosition.x + 75 &&
			sf::Mouse::getPosition(window).y > loadModelButtonPosition.y - 20 && sf::Mouse::getPosition(window).y < loadModelButtonPosition.y + 20) {
			loadModelButton.setFillColor(sf::Color(100, 100, 100));

			//if clicked, load the model
			if (leftMouseClickedThisFrame) {
				sf::Image img = LoadImageThroughDialog();

				//check if img is valid
				if (img.getSize().x > 0) {
					Obstacle obs(img, 0.25, 0.5, 0.7);
					
					fluidSim.simParams.obstacles.clear();
					fluidSim.simParams.obstacles.push_back(obs);


					fluidSim.UpdateSField();
					fluidSim.SetupWindTunnelBoundaries();
					
				}
				

			}

			

			
			
		}
		else {
			loadModelButton.setFillColor(sf::Color(70, 70, 70));
		}


		//if the mouse is over the button, change the color
		if (sf::Mouse::getPosition(window).x > resetButtonPosition.x - 75 && sf::Mouse::getPosition(window).x < resetButtonPosition.x + 75 &&
			sf::Mouse::getPosition(window).y > resetButtonPosition.y - 20 && sf::Mouse::getPosition(window).y < resetButtonPosition.y + 20) {
			resetButton.setFillColor(sf::Color(100, 100, 100));

			//if clicked, reset the simulation
			if (leftMouseClickedThisFrame) {
				fluidSim.Reset();
				//fluidSim.SetupWindTunnelBoundaries(4);

				fluidSim.UpdateSField();
				fluidSim.SetupWindTunnelBoundaries();


				//fluidSim.SetObstacle(0.25, 0.5, 0.08);

				//fluidSim.SetObstacle(0.60, 0.5, 0.08);
			}
		}
		else {
			resetButton.setFillColor(sf::Color(70, 70, 70));
		}

		
		//save the currnet image to up 2 folders and into "images"
		//image.saveToFile("C:\\Users\\aspen\\Desktop\\Euler-Fluid-Simulation\\Images\\image" + std::to_string(displayParams.frameCount) + ".png");

		displayParams.frameCount += 1;


		//draw the window
		window.display();
		



	}



}
