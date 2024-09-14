#include <SFML/Graphics.hpp>
#include <windows.h>
#include <commdlg.h>
#include <string>
#include <iostream>

// Helper function to convert UTF-16 (wchar_t) to UTF-8 (std::string)
std::string utf16ToUtf8(const std::wstring& wstr)
{
    if (wstr.empty()) return std::string();
    int size_needed = WideCharToMultiByte(
        CP_UTF8, 0, &wstr[0], (int)wstr.size(),
        NULL, 0, NULL, NULL
    );
    std::string strTo(size_needed, 0);
    WideCharToMultiByte(
        CP_UTF8, 0, &wstr[0], (int)wstr.size(),
        &strTo[0], size_needed, NULL, NULL
    );
    return strTo;
}

sf::Image LoadImageThroughDialog()
{
    sf::Image image;

    OPENFILENAMEW ofn;          // Common dialog box structure
    WCHAR szFile[MAX_PATH] = { 0 };     // Buffer for file name, initialized to zero

    // Initialize OPENFILENAME
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;       // If you have a window handle, put it here
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = MAX_PATH;

    // Set the file filters
    ofn.lpstrFilter = L"Image Files\0*.bmp;*.jpg;*.jpeg;*.png;*.tga\0All Files\0*.*\0";
    ofn.nFilterIndex = 1;

    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    // Display the Open dialog box
    if (GetOpenFileNameW(&ofn) == TRUE)
    {
        // Convert the wide string to UTF-8
        std::wstring ws(ofn.lpstrFile);
        std::string filePath = utf16ToUtf8(ws);

        // Load the image from the selected file
        if (!image.loadFromFile(filePath))
        {
            // Handle loading error
            MessageBoxW(NULL, L"Failed to load image.", L"Error", MB_OK | MB_ICONERROR);
        }
    }
    else
    {
        // Get error from the common dialog box
        DWORD error = CommDlgExtendedError();

        // Print or display more detailed error information
        if (error == 0)
        {
            //MessageBoxW(NULL, L"User canceled the dialog.", L"Information", MB_OK | MB_ICONINFORMATION);
			return image;
        }
        else
        {
            // Print out the error code for further investigation
            std::wstring errorMessage = L"An error occurred. Error code: " + std::to_wstring(error);
            MessageBoxW(NULL, errorMessage.c_str(), L"Error", MB_OK | MB_ICONERROR);
        }
    }

    return image;
}
