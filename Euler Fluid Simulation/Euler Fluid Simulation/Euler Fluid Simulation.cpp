#include <iostream>
#include <vector>
#include "FluidSim.hpp"
#include <thread>
#include "SFML/Graphics.hpp"
#include <string>
#include <chrono>
#include "Helvetica.hpp"
#include "GUI.hpp"

int main() {
    SimParameters simParams;
    DisplayParameters displayParams;

    size_t fixedHeight = 800;
    simParams.SetGridSize(300, 200);
    simParams.numIterations = 200;

#ifdef CXXDROID_COMPAT
    simParams.SetGridSize(220, 100);
    simParams.numIterations = 100;
    fixedHeight = 1080;
#endif

    displayParams.windowWidth = (size_t)(fixedHeight * (float)simParams.n_x / (float)simParams.n_y);
    displayParams.windowHeight = fixedHeight;
    displayParams.maxFps = 60;

    simParams.windTunnelSpeed = 2.f;
    simParams.overRelaxation = 1.9f;
    simParams.gravity = 10.f;

	sf::RenderWindow window(sf::VideoMode(displayParams.windowWidth, displayParams.windowHeight), "Euler Fluid Simulation", sf::Style::Close);
    window.setFramerateLimit(displayParams.maxFps);

    debugDrawWindow = &window;
    FluidSim fluidSim(simParams, displayParams);
    //Obstacle obs(0.30, 0.48, 0.13, 0.13);

    //fluidSim.displayParams.displayMode = DisplayParameters::DisplayMode::DYE;
    //fluidSim.displayParams.displayMode = DisplayParameters::DisplayMode::PRESSURE;
    //fluidSim.displayParams.displayMode = DisplayParameters::DisplayMode::VELOCITIES;
    fluidSim.displayParams.displayMode = DisplayParameters::DisplayMode::VORTICITY;


    // Adding dye sources
    float dyeLineSpacing = 1.f / (float)simParams.n_y;
    for (size_t i = 0; i < simParams.n_y; i++) {
        if (i % 10 == 0) {
            fluidSim.AddDyeSource(0.01f, 0.05f + i * dyeLineSpacing, 0.02f, 0.001f, 0.f);
        }
    }




    GUI gui(&fluidSim, displayParams.windowWidth, displayParams.windowHeight, displayParams.maxFps);


    gui.LoadBuiltInModel(0);

    sf::Clock frameClock, simClock, renderClock;

    while (window.isOpen()) {
        bool leftMouseClickedThisFrame = false;
        bool mouseUnclickedThisFrame = false;
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
                window.close();
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left)
                leftMouseClickedThisFrame = true;

            if (event.type ==sf::Event::TouchBegan)
				leftMouseClickedThisFrame = true;

            if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left)
				mouseUnclickedThisFrame = true;
            if (event.type == sf::Event::TouchEnded)
                mouseUnclickedThisFrame = true;
        }

        window.clear();

#ifndef CXXDROID_COMPAT
        sf::Vector2i pointerPositionInt = sf::Mouse::getPosition(window);
#endif
#ifdef CXXDROID_COMPAT
        sf::Vector2i pointerPositionInt = sf::Touch::getPosition(0);
#endif



        sf::Vector2f pointerPosition = sf::Vector2f((float)pointerPositionInt.x, (float)pointerPositionInt.y);

        float currentTime = frameClock.restart().asSeconds();
        float fps = 1.f / (currentTime);

        simClock.restart();
        fluidSim.Simulate(0.005f);
        float simTime = (float)simClock.restart().asMicroseconds();

        renderClock.restart();
        fluidSim.Render(window);
        float renderTime = (float)renderClock.restart().asMicroseconds();

        gui.updateText(fps, simTime, renderTime, currentTime);
        gui.update(window, pointerPosition, leftMouseClickedThisFrame, mouseUnclickedThisFrame);
        gui.draw(window);

        window.display();
    }
}
