#include "GUI.hpp"
#include <iostream>
#include "Helvetica.hpp"
#include "Airfoil_AA_30_Degrees.hpp"

#ifndef CXXDROID_COMPAT
#include "ModelLoading.hpp"
#endif

GUI::GUI(FluidSim* fluidSim, int windowWidth, int windowHeight, int maxFps) {
    this->fluidSim = fluidSim;
    selectedModelIndex = 0;  // Start with the first model

    // Load font from memory
    if (!font.loadFromMemory(helveticaFontBytes, 714756)) {
        std::cerr << "Error loading font" << std::endl;
    }

    // Initialize texts
    text.setFont(font);
    text.setCharacterSize(18);
    text.setFillColor(sf::Color::White);
    text.setPosition(5, 5);

    text2.setFont(font);
    text2.setCharacterSize(18);
    text2.setFillColor(sf::Color::White);
    text2.setPosition(320, 5);

    text3.setFont(font);
    text3.setCharacterSize(22);
    text3.setFillColor(sf::Color::White);
    text3.setPosition(1200, 5);

	copyrightText.setFont(font);
	copyrightText.setCharacterSize(28);
    copyrightText.setFillColor(sf::Color(255, 255, 255, 140));
    copyrightText.setString(sf::String(L"©Aspen Erlandsson 2024"));
	copyrightText.setPosition(5, (float)(windowHeight - 34));
    

    // Initialize buttons
    loadModelButtonPosition = sf::Vector2f((float)(windowWidth - 100), 50.f);
    loadModelButton.setSize(sf::Vector2f(150, 40));
    loadModelButton.setFillColor(sf::Color(70, 70, 70));
    loadModelButton.setOrigin(loadModelButton.getLocalBounds().width / 2, loadModelButton.getLocalBounds().height / 2);
    loadModelButton.setPosition(loadModelButtonPosition);

    loadModelText.setFont(font);
    loadModelText.setCharacterSize(24);
    loadModelText.setFillColor(sf::Color::White);
    loadModelText.setString("Load PNG");
    loadModelText.setOrigin(loadModelText.getLocalBounds().width / 2, loadModelText.getLocalBounds().height / 2 + 5);
    loadModelText.setPosition(loadModelButtonPosition);

    resetButtonPosition = sf::Vector2f((float)(windowWidth - 100), 100.f);
    resetButton.setSize(sf::Vector2f(150, 40));
    resetButton.setFillColor(sf::Color(70, 70, 70));
    resetButton.setOrigin(resetButton.getLocalBounds().width / 2, resetButton.getLocalBounds().height / 2);
    resetButton.setPosition(resetButtonPosition);

    resetText.setFont(font);
    resetText.setCharacterSize(24);
    resetText.setFillColor(sf::Color::White);
    resetText.setString("Reset Dye");
    resetText.setOrigin(resetText.getLocalBounds().width / 2, resetText.getLocalBounds().height / 2 + 5);
    resetText.setPosition(resetButtonPosition);

    // Initialize model title and buttons
    modelTitleText.setFont(font);
    modelTitleText.setCharacterSize(18);
    modelTitleText.setFillColor(sf::Color::White);
    modelTitleText.setString("Built-in Models:");
    modelTitleText.setPosition((float)(windowWidth - 160), (float)150);  // Title position

    selectedModelText.setFont(font);
    selectedModelText.setCharacterSize(18);
    selectedModelText.setFillColor(sf::Color(245, 81, 81));
    selectedModelText.setString(modelNames[selectedModelIndex]);
    selectedModelText.setPosition((float)(windowWidth - 100), (float)178);  // Selected model position
    //center
    selectedModelText.setOrigin(selectedModelText.getLocalBounds().width / 2, selectedModelText.getLocalBounds().height / 2);


    previousModelButtonPosition = sf::Vector2f((float)(windowWidth - 180), 200.f);
    previousModelButton.setSize(sf::Vector2f(70, 50));
    previousModelButton.setFillColor(sf::Color(70, 70, 70));
    previousModelButton.setPosition(previousModelButtonPosition);

    previousModelText.setFont(font);
    previousModelText.setCharacterSize(24);
    previousModelText.setFillColor(sf::Color::White);
    previousModelText.setString("<");
    previousModelText.setOrigin(previousModelText.getLocalBounds().width / 2, previousModelText.getLocalBounds().height / 2);
    previousModelText.setPosition(previousModelButtonPosition.x + 32, previousModelButtonPosition.y + 12);

    nextModelButtonPosition = sf::Vector2f((float)(windowWidth - 90), 200.f);
    nextModelButton.setSize(sf::Vector2f(70, 50));
    nextModelButton.setFillColor(sf::Color(70, 70, 70));
    nextModelButton.setPosition(nextModelButtonPosition);

    nextModelText.setFont(font);
    nextModelText.setCharacterSize(24);
    nextModelText.setFillColor(sf::Color::White);
    nextModelText.setString(">");
    nextModelText.setOrigin(nextModelText.getLocalBounds().width / 2, nextModelText.getLocalBounds().height / 2);
    nextModelText.setPosition(nextModelButtonPosition.x + 32, nextModelButtonPosition.y + 12);




    // Initialize visualization title and buttons
    visualizationTitleText.setFont(font);
    visualizationTitleText.setCharacterSize(18);
    visualizationTitleText.setFillColor(sf::Color::White);
    visualizationTitleText.setString("Visualization:");
    visualizationTitleText.setPosition((float)(windowWidth - 160), (float)300);  // Title position

    selectedVisualizationText.setFont(font);
    selectedVisualizationText.setCharacterSize(18);
    selectedVisualizationText.setFillColor(sf::Color(245, 81, 81));
    selectedVisualizationText.setString(visualizationNames[selectedVisualizationIndex]);
    selectedVisualizationText.setPosition((float)(windowWidth - 100), (float)328);  // Selected visualization position
    selectedVisualizationText.setOrigin(selectedVisualizationText.getLocalBounds().width / 2, selectedVisualizationText.getLocalBounds().height / 2);

    previousVisualizationButton.setSize(sf::Vector2f(70, 50));
    previousVisualizationButton.setFillColor(sf::Color(70, 70, 70));
    previousVisualizationButton.setPosition((float)(windowWidth - 180), 350.f);

    previousVisualizationText.setFont(font);
    previousVisualizationText.setCharacterSize(24);
    previousVisualizationText.setFillColor(sf::Color::White);
    previousVisualizationText.setString("<");
    previousVisualizationText.setOrigin(previousVisualizationText.getLocalBounds().width / 2, previousVisualizationText.getLocalBounds().height / 2);
    previousVisualizationText.setPosition(previousVisualizationButton.getPosition().x + 32, previousVisualizationButton.getPosition().y + 12);

    nextVisualizationButton.setSize(sf::Vector2f(70, 50));
    nextVisualizationButton.setFillColor(sf::Color(70, 70, 70));
    nextVisualizationButton.setPosition((float)(windowWidth - 90), 350.f);

    nextVisualizationText.setFont(font);
    nextVisualizationText.setCharacterSize(24);
    nextVisualizationText.setFillColor(sf::Color::White);
    nextVisualizationText.setString(">");
    nextVisualizationText.setOrigin(nextVisualizationText.getLocalBounds().width / 2, nextVisualizationText.getLocalBounds().height / 2);
    nextVisualizationText.setPosition(nextVisualizationButton.getPosition().x + 32, nextVisualizationButton.getPosition().y + 12);

    
}

void GUI::update(sf::RenderWindow& window, const sf::Vector2f& pointerPosition, bool leftMouseClickedThisFrame, bool mouseUnclickedThisFrame) {
    // Handle load model button interaction (same as before)
    if (pointerPosition.x > loadModelButtonPosition.x - 75 && pointerPosition.x < loadModelButtonPosition.x + 75 &&
        pointerPosition.y > loadModelButtonPosition.y - 20 && pointerPosition.y < loadModelButtonPosition.y + 20) {
        loadModelButton.setFillColor(sf::Color(100, 100, 100));
        if (leftMouseClickedThisFrame) {
#ifndef CXXDROID_COMPAT
            sf::Image img = LoadImageThroughDialog();

#endif
#ifdef CXXDROID_COMPAT
            sf::Image img;
#endif
            if (img.getSize().x > 0) {
                Obstacle obs(img, 0.7f, 0.5f, 0.7f);
                fluidSim->simParams.obstacles.clear();
                fluidSim->simParams.obstacles.push_back(obs);
                fluidSim->UpdateSField();
                fluidSim->SetupWindTunnelBoundaries();
            }
        }
    }
    else {
        loadModelButton.setFillColor(sf::Color(70, 70, 70));
    }

    // Handle reset button interaction (same as before)
    if (pointerPosition.x > resetButtonPosition.x - 75 && pointerPosition.x < resetButtonPosition.x + 75 &&
        pointerPosition.y > resetButtonPosition.y - 20 && pointerPosition.y < resetButtonPosition.y + 20) {
        resetButton.setFillColor(sf::Color(100, 100, 100));
        if (leftMouseClickedThisFrame) {
            std::fill(fluidSim->mField, fluidSim->mField + fluidSim->simParams.n_cells, 1.0f);
            fluidSim->UpdateSField();
            fluidSim->SetupWindTunnelBoundaries();
        }
    }
    else {
        resetButton.setFillColor(sf::Color(70, 70, 70));
    }

    // Handle previous model button interaction
    if (pointerPosition.x > previousModelButtonPosition.x && pointerPosition.x < previousModelButtonPosition.x + 70 &&
        pointerPosition.y > previousModelButtonPosition.y && pointerPosition.y < previousModelButtonPosition.y + 50) {
        previousModelButton.setFillColor(sf::Color(100, 100, 100));
        if (leftMouseClickedThisFrame) {
            selectedModelIndex = (selectedModelIndex - 1 + 3) % 3;  // Cycle backwards
            std::cout << "Selected model: " << modelNames[selectedModelIndex] << std::endl;
            LoadBuiltInModel(selectedModelIndex);
            selectedModelText.setString(modelNames[selectedModelIndex]);
            selectedModelText.setOrigin(selectedModelText.getLocalBounds().width / 2, selectedModelText.getLocalBounds().height / 2);

        }
    }
    else {
        previousModelButton.setFillColor(sf::Color(70, 70, 70));
    }

    // Handle next model button interaction
    if (pointerPosition.x > nextModelButtonPosition.x && pointerPosition.x < nextModelButtonPosition.x + 70 &&
        pointerPosition.y > nextModelButtonPosition.y && pointerPosition.y < nextModelButtonPosition.y + 50) {
        nextModelButton.setFillColor(sf::Color(100, 100, 100));
        if (leftMouseClickedThisFrame) {
            selectedModelIndex = (selectedModelIndex + 1) % 3;  // Cycle forwards
            std::cout << "Selected model: " << modelNames[selectedModelIndex] << std::endl;
            LoadBuiltInModel(selectedModelIndex);
            selectedModelText.setString( modelNames[selectedModelIndex]);
            selectedModelText.setOrigin(selectedModelText.getLocalBounds().width / 2, selectedModelText.getLocalBounds().height / 2);
        }
    }
    else {
        nextModelButton.setFillColor(sf::Color(70, 70, 70));
    }


    // Handle previous visualization button interaction
    if (pointerPosition.x > previousVisualizationButton.getPosition().x && pointerPosition.x < previousVisualizationButton.getPosition().x + 70 &&
        pointerPosition.y > previousVisualizationButton.getPosition().y && pointerPosition.y < previousVisualizationButton.getPosition().y + 50) {
        previousVisualizationButton.setFillColor(sf::Color(100, 100, 100));
        if (leftMouseClickedThisFrame) {
            selectedVisualizationIndex = (selectedVisualizationIndex - 1 + 4) % 4;  // Cycle backwards
            std::cout << "Selected visualization: " << visualizationNames[selectedVisualizationIndex] << std::endl;
            fluidSim->displayParams.displayMode = static_cast<DisplayParameters::DisplayMode>(selectedVisualizationIndex);
            selectedVisualizationText.setString(visualizationNames[selectedVisualizationIndex]);
            selectedVisualizationText.setOrigin(selectedVisualizationText.getLocalBounds().width / 2, selectedVisualizationText.getLocalBounds().height / 2);
        }
    }
    else {
        previousVisualizationButton.setFillColor(sf::Color(70, 70, 70));
    }

    // Handle next visualization button interaction
    if (pointerPosition.x > nextVisualizationButton.getPosition().x && pointerPosition.x < nextVisualizationButton.getPosition().x + 70 &&
        pointerPosition.y > nextVisualizationButton.getPosition().y && pointerPosition.y < nextVisualizationButton.getPosition().y + 50) {
        nextVisualizationButton.setFillColor(sf::Color(100, 100, 100));
        if (leftMouseClickedThisFrame) {
            selectedVisualizationIndex = (selectedVisualizationIndex + 1) % 4;  // Cycle forwards
            std::cout << "Selected visualization: " << visualizationNames[selectedVisualizationIndex] << std::endl;
            fluidSim->displayParams.displayMode = static_cast<DisplayParameters::DisplayMode>(selectedVisualizationIndex);
            selectedVisualizationText.setString(visualizationNames[selectedVisualizationIndex]);
            selectedVisualizationText.setOrigin(selectedVisualizationText.getLocalBounds().width / 2, selectedVisualizationText.getLocalBounds().height / 2);
        }
    }
    else {
        nextVisualizationButton.setFillColor(sf::Color(70, 70, 70));
    }



    sf::Vector2f simSpaceMousePosition = sf::Vector2f((float)pointerPosition.x / window.getSize().y, (float)pointerPosition.y / window.getSize().y);

    sf::Vector2f mouseVelocity = simSpaceMousePosition - lastMousePosition;

    float lastMouseUpdateTime = mouseSpeedClock.getElapsedTime().asSeconds();

    mouseVelocity = mouseVelocity / lastMouseUpdateTime * 0.0f;





    bool mouseOverBoundary = fluidSim->PositionIsBoundary(simSpaceMousePosition.x, simSpaceMousePosition.y);


    if (leftMouseClickedThisFrame && mouseOverBoundary)
    {
        isDragging = true;
        dragStartOffset = simSpaceMousePosition - sf::Vector2f(fluidSim->simParams.obstacles[0].x, fluidSim->simParams.obstacles[0].y);

    }

    if (mouseUnclickedThisFrame)
    {
        isDragging = false;

    }


    if (isDragging)
    {
        fluidSim->simParams.mouseVelocity = mouseVelocity;
    }
    else
    {
		fluidSim->simParams.mouseVelocity = sf::Vector2f(0, 0);
	}

    //if dragging, update the obstacle position and initiate refresh of the simulation
    if (isDragging)
    {
		fluidSim->simParams.obstacles[0].x = simSpaceMousePosition.x - dragStartOffset.x;
		fluidSim->simParams.obstacles[0].y = simSpaceMousePosition.y - dragStartOffset.y;
		fluidSim->UpdateSField();
		fluidSim->SetupWindTunnelBoundaries();
       
	}

    lastMousePosition = pointerPosition;
    mouseSpeedClock.restart();


}

void GUI::draw(sf::RenderWindow& window) {
    // Draw elements
    window.draw(text);
    window.draw(text2);
    window.draw(text3);
    window.draw(loadModelButton);
    window.draw(loadModelText);
    window.draw(resetButton);
    window.draw(resetText);
    window.draw( selectedModelText);
    window.draw(copyrightText);

    // Draw built-in models title and buttons
    window.draw(modelTitleText);
    window.draw(previousModelButton);
    window.draw(previousModelText);
    window.draw(nextModelButton);
    window.draw(nextModelText);

    // Draw visualization selector
    window.draw(visualizationTitleText);
    window.draw(previousVisualizationButton);
    window.draw(previousVisualizationText);
    window.draw(nextVisualizationButton);
    window.draw(nextVisualizationText);
    window.draw(selectedVisualizationText);

}



void GUI::updateText(float fps, float simTime, float renderTime, float currentTime) {
    // Update the text strings
    text.setString("GRID: " + std::to_string(fluidSim->simParams.n_x) + "x" + std::to_string(fluidSim->simParams.n_y) + ", FPS: " + std::to_string(fps));
    text2.setString("Fluid Sim Time: " + std::to_string(simTime / 1000.f) + "ms, Fluid Render Time: " + std::to_string(renderTime / 1000.f) + "ms" + ", Frame Time: " + std::to_string(currentTime * 1000) + "ms");
    text3.setString("Avg Divergence: " + std::to_string(fluidSim->simStats.avgVelocityFieldDivergence) + ", Max Divergence: " + std::to_string(fluidSim->simStats.maxVelocityFieldDivergence));
}


void GUI::LoadBuiltInModel(int index)
{
    //if name is circle
    if (modelNames[index] == "Circle")
    {
		Obstacle obs(0.5f, 0.5f, 0.14f);
		fluidSim->simParams.obstacles.clear();
		fluidSim->simParams.obstacles.push_back(obs);
		fluidSim->UpdateSField();
		fluidSim->SetupWindTunnelBoundaries();
	}

    //if name is square
    if (modelNames[index] == "Square")
    {
		Obstacle obs(0.5f, 0.5f, 0.2f, 0.2f);
		fluidSim->simParams.obstacles.clear();
		fluidSim->simParams.obstacles.push_back(obs);
		fluidSim->UpdateSField();
		fluidSim->SetupWindTunnelBoundaries();
	}

	//if name is airfoil
    if (modelNames[index] == "Airfoil High AA")
    {

        loadedImage.loadFromMemory(airfoil_aa_30_degrees, 118227);

        Obstacle obs(loadedImage, 0.8f, 0.5f, 0.8f);
        
        fluidSim->simParams.obstacles.clear();
		fluidSim->simParams.obstacles.push_back(obs);
		fluidSim->UpdateSField();
		fluidSim->SetupWindTunnelBoundaries();
	}




}