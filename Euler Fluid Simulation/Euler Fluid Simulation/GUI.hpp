#pragma once

#include "SFML/Graphics.hpp"
#include "FluidSim.hpp"

class GUI {
public:
    // Constructor
    GUI(FluidSim* fluidSim, int windowWidth, int windowHeight, int maxFps);

    // Methods
    void update(sf::RenderWindow& window, const sf::Vector2f& pointerPosition, bool leftMouseClickedThisFrame, bool mouseUnclickedThisFrame);
    void draw(sf::RenderWindow& window);
    void updateText(float fps, float simTime, float renderTime, float currentTime);

    void LoadBuiltInModel(int index);

private:
    FluidSim* fluidSim;  // Pointer to FluidSim object

    // GUI elements
    sf::Text text, text2, text3;
    sf::RectangleShape loadModelButton, resetButton;
    sf::Text loadModelText, resetText;

    sf::RectangleShape previousModelButton, nextModelButton;
    sf::Text previousModelText, nextModelText;
    sf::Text modelTitleText;

    sf::Text copyrightText;

    sf::Text selectedModelText;

    sf::Font font;

    sf::Vector2f dragBeginPosition;

    sf::Image loadedImage;

    // Button positions
    sf::Vector2f loadModelButtonPosition;
    sf::Vector2f resetButtonPosition;
    sf::Vector2f previousModelButtonPosition;
    sf::Vector2f nextModelButtonPosition;

    sf::Vector2f dragStartOffset;


    sf::Vector2f lastMousePosition;

    sf::Clock mouseSpeedClock;

    bool isDragging = false;

    // Built-in model tracking
    int selectedModelIndex;
    std::string modelNames[3] = { "Circle", "Square", "Airfoil High AA" };
};
