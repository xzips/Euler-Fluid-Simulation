// Credit to Kenneth Moreland: https://www.kennethmoreland.com/color-advice/ for the color map values
//(implementation by me still, but specific values were taken from Kenneth Moreland's website)

#pragma once


#include <SFML/Graphics/Color.hpp>
#include <vector>
#include <algorithm>


namespace SciColorMaps
{

    struct ColorPoint {
        float scalar;
        sf::Color color;
    };

    sf::Color InterpolateColor(const sf::Color& c1, const sf::Color& c2, float t) {
        return sf::Color(
            static_cast<sf::Uint8>(c1.r + t * (c2.r - c1.r)),
            static_cast<sf::Uint8>(c1.g + t * (c2.g - c1.g)),
            static_cast<sf::Uint8>(c1.b + t * (c2.b - c1.b))
        );
    }

    sf::Color GetColorFromMap(float s, const std::vector<ColorPoint>& colorMap) {
        // Clamp the scalar to the range [0, 1]
        s = std::clamp(s, 0.0f, 1.0f);

        // Find two adjacent points in the color map to interpolate between
        for (size_t i = 0; i < colorMap.size() - 1; ++i) {
            if (s >= colorMap[i].scalar && s <= colorMap[i + 1].scalar) {
                // Interpolation factor
                float t = (s - colorMap[i].scalar) / (colorMap[i + 1].scalar - colorMap[i].scalar);
                return InterpolateColor(colorMap[i].color, colorMap[i + 1].color, t);
            }
        }

        return colorMap.back().color;
    }

    static const std::vector<ColorPoint> GetViridisColorMap() {
        return {
            {0.0f, sf::Color(68, 1, 84)},
            {0.06666666667f, sf::Color(72, 26, 108)},
            {0.13333333333f, sf::Color(71, 47, 125)},
            {0.2f, sf::Color(65, 68, 135)},
            {0.26666666667f, sf::Color(57, 86, 140)},
            {0.33333333333f, sf::Color(49, 104, 142)},
            {0.4f, sf::Color(42, 120, 142)},
            {0.46666666667f, sf::Color(35, 136, 142)},
            {0.53333333333f, sf::Color(31, 152, 139)},
            {0.6f, sf::Color(34, 168, 132)},
            {0.66666666667f, sf::Color(53, 183, 121)},
            {0.73333333333f, sf::Color(84, 197, 104)},
            {0.8f, sf::Color(122, 209, 81)},
            {0.86666666667f, sf::Color(165, 219, 54)},
            {0.93333333333f, sf::Color(210, 226, 27)},
            {1.0f, sf::Color(253, 231, 37)}
        };
    }


    // New Plasma Color Map
    static const std::vector<ColorPoint> GetPlasmaColorMap() {
        return {
            {0.0f, sf::Color(13, 8, 135)},
            {0.06666666667f, sf::Color(51, 5, 151)},
            {0.13333333333f, sf::Color(80, 2, 162)},
            {0.2f, sf::Color(106, 0, 168)},
            {0.26666666667f, sf::Color(132, 5, 167)},
            {0.33333333333f, sf::Color(156, 23, 158)},
            {0.4f, sf::Color(177, 42, 144)},
            {0.46666666667f, sf::Color(195, 61, 128)},
            {0.53333333333f, sf::Color(211, 81, 113)},
            {0.6f, sf::Color(225, 100, 98)},
            {0.66666666667f, sf::Color(237, 121, 83)},
            {0.73333333333f, sf::Color(246, 143, 68)},
            {0.8f, sf::Color(252, 166, 54)},
            {0.86666666667f, sf::Color(254, 192, 41)},
            {0.93333333333f, sf::Color(249, 220, 36)},
            {1.0f, sf::Color(240, 249, 33)}
        };
    }

    // New Blackbody Color Map
    static const std::vector<ColorPoint> GetBlackbodyColorMap() {
        return {
            {0.0f, sf::Color(0, 0, 0)},
            {0.06666666667f, sf::Color(36, 15, 9)},
            {0.13333333333f, sf::Color(62, 22, 17)},
            {0.2f, sf::Color(90, 27, 22)},
            {0.26666666667f, sf::Color(119, 30, 26)},
            {0.33333333333f, sf::Color(150, 33, 30)},
            {0.4f, sf::Color(180, 38, 34)},
            {0.46666666667f, sf::Color(197, 65, 28)},
            {0.53333333333f, sf::Color(214, 88, 19)},
            {0.6f, sf::Color(228, 112, 7)},
            {0.66666666667f, sf::Color(231, 141, 18)},
            {0.73333333333f, sf::Color(233, 169, 29)},
            {0.8f, sf::Color(233, 195, 39)},
            {0.86666666667f, sf::Color(231, 222, 50)},
            {0.93333333333f, sf::Color(246, 240, 144)},
            {1.0f, sf::Color(255, 255, 255)}
        };
    }
    // New Inferno Color Map
    static const std::vector<ColorPoint> GetInfernoColorMap() {
        return {
            {0.0f, sf::Color(0, 0, 4)},
            {0.06666666667f, sf::Color(12, 8, 38)},
            {0.13333333333f, sf::Color(36, 12, 79)},
            {0.2f, sf::Color(66, 10, 104)},
            {0.26666666667f, sf::Color(93, 18, 110)},
            {0.33333333333f, sf::Color(120, 28, 109)},
            {0.4f, sf::Color(147, 38, 103)},
            {0.46666666667f, sf::Color(174, 48, 92)},
            {0.53333333333f, sf::Color(199, 62, 76)},
            {0.6f, sf::Color(221, 81, 58)},
            {0.66666666667f, sf::Color(237, 105, 37)},
            {0.73333333333f, sf::Color(248, 133, 15)},
            {0.8f, sf::Color(252, 165, 10)},
            {0.86666666667f, sf::Color(250, 198, 45)},
            {0.93333333333f, sf::Color(242, 230, 97)},
            {1.0f, sf::Color(252, 255, 164)}
        };
    }


    static const std::vector<ColorPoint> GetFastColorMap() {
        return {
            {0.0f, sf::Color(14, 14, 120)},
            {0.06666666667f, sf::Color(39, 55, 153)},
            {0.13333333333f, sf::Color(55, 94, 187)},
            {0.2f, sf::Color(70, 133, 215)},
            {0.26666666667f, sf::Color(85, 171, 234)},
            {0.33333333333f, sf::Color(118, 202, 241)},
            {0.4f, sf::Color(159, 226, 237)},
            {0.46666666667f, sf::Color(205, 239, 215)},
            {0.53333333333f, sf::Color(236, 230, 171)},
            {0.6f, sf::Color(244, 207, 125)},
            {0.66666666667f, sf::Color(240, 177, 97)},
            {0.73333333333f, sf::Color(231, 146, 73)},
            {0.8f, sf::Color(215, 114, 54)},
            {0.86666666667f, sf::Color(197, 82, 40)},
            {0.93333333333f, sf::Color(174, 54, 35)},
            {1.0f, sf::Color(150, 20, 30)}
        };


    }


    //"Viridis is a perceptually uniform color map with monotonically increasing luminance and a pleasant smooth arc through blue, green, and yellow hues.
    sf::Color Viridis(float s) {
        static const std::vector<ColorPoint> viridisMap = GetViridisColorMap();
        return GetColorFromMap(s, viridisMap);
    }

    //"Plasma is a perceptually uniform color map with monotonically increasing luminance and a pleasant smooth arc through blue, purple, and yellow hues."
    sf::Color Plasma(float s) {
        static const std::vector<ColorPoint> plasmaMap = GetPlasmaColorMap();
        return GetColorFromMap(s, plasmaMap);
    }


    //The black body color map is based on colors from black-body radiation. The colors are are not exact to those of black-body radiation but are designed to have a constant increase in brightness throughout.
    sf::Color Blackbody(float s) {
        static const std::vector<ColorPoint> blackbodyMap = GetBlackbodyColorMap();
        return GetColorFromMap(s, blackbodyMap);
    }


    // "Inferno is a high-contrast, perceptually uniform color map with a smooth transition from dark blue to bright yellow, ideal for scientific visualization."
    sf::Color Inferno(float s) {
        static const std::vector<ColorPoint> infernoMap = GetInfernoColorMap();
        return GetColorFromMap(s, infernoMap);
    }

    //This color map was designed to be the default colormap in the ParaView scientific visualization tool. It is a diverging (double-ended) color map with a smooth transition in the middle to prevent artifacts at the midpoint. This colormap is designed to be used on 3D surfaces, so it avoids getting too dark at the ends.
    sf::Color Fast(float s) {
        static const std::vector<ColorPoint> fastMap = GetFastColorMap();
        return GetColorFromMap(s, fastMap);
    }
}