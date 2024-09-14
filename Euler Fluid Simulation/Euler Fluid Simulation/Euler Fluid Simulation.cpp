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
	sf::RenderWindow window(sf::VideoMode(800, 800), "Euler Fluid Simulation");

	//framerate limit to 60
	window.setFramerateLimit(30);

	//create RenderTexture
	sf::RenderTexture renderTexture;
	renderTexture.create(800, 800);



	//fluidSim.SetObstacle(16, 16, 10);

	
	
	//create text object for drawing cell values
	
	sf::Font font;
	if (!font.loadFromFile("C:\\Users\\aspen\\Desktop\\SFMLFonts\\Helvetica.ttf"))
	{
		std::cout << "Error loading font" << std::endl;
	}

	sf::Text text;
	text.setFont(font);
	text.setCharacterSize(10);

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




		float inVel = 2;
		
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



		fluidSim.SetObstacle(0.35, 0.5, 0.2);




		float pipeH = 0.1 * simParams.n_y;
		size_t minJ = std::floor(0.5 * simParams.n_y - 0.5 * pipeH);
		size_t maxJ = std::floor(0.5 * simParams.n_y + 0.5 * pipeH);

		for (size_t j = minJ; j < maxJ; j++)
			fluidSim.mField[2 * simParams.n_y + j] = 0.0;







		
		float min_P = INFINITY;
		float max_P = -INFINITY;

		//find min and max
		for (size_t i = 1; i < simParams.n_x; i++)
		{
			for (size_t j = 1; j < simParams.n_y; j++)
			{
				if (fluidSim.mField[i * simParams.n_y + j] < min_P)
					min_P = fluidSim.mField[i * simParams.n_y + j];
				if (fluidSim.mField[i * simParams.n_y + j] > max_P)
					max_P = fluidSim.mField[i * simParams.n_y + j];
			}
		}

		std::cout << "min P: " << min_P << " max P: " << max_P << std::endl;

		//sample V field and draw it
		for (size_t i = 0; i < simParams.n_x; i++)
		{
			for (size_t j = 0; j < simParams.n_y; j++)
			{

				float P_norm = 1- (fluidSim.mField[i * simParams.n_y + j] - min_P) / (max_P - min_P);

				//draw the pressure field as a rectangle with the color of the pressure value
				sf::RectangleShape rect(sf::Vector2f(800.f / simParams.n_x, 800.f / simParams.n_y));
				rect.setPosition(i * 800.f / simParams.n_x, j * 800.f / simParams.n_y);
				rect.setFillColor(sf::Color(255 * P_norm, 255 * P_norm, 255 * P_norm, 255));

				//if S field greater than 0 change rect color to blue
				if (fluidSim.sField[i * simParams.n_y + j] == 0)
					rect.setFillColor(sf::Color(0, 0, 255, 255));



				window.draw(rect);

				//set text position to the cell
				//text.setPosition(i * 800.f / simParams.n_x, j * 800.f / simParams.n_y);
				//set the text string to the value of the cell
				//std::string dispStr = std::to_string((long int)fluidSim.SampleField(i, j, FieldType::V_FIELD));

				//p field
				//std::string dispStr = std::to_string(P_norm);
				
				//text.setString(dispStr);
				//draw the text
				//window.draw(text);


			}
		}
		

	

		//simulate the fluid
		fluidSim.Simulate(0.01f);




		//draw the window
		window.display();
		



	}



}
