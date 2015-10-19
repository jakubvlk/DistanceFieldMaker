//
//  main.cpp
//  DistanceFieldMaker
//
//  Created by Jakub Vlk on 21/09/15.
//  Copyright (c) 2015 Jakub Vlk. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <float.h>
#include <math.h>
#include <ctime>
#include <array>
#include "lz4.h"



using namespace::std;
typedef unsigned char byte;

static int version;
static int depth, height, width;
static int size;

vector<byte> voxelsVec;
vector<byte> originalVoxelsVec;
vector<double> distanceFieldVec;

static float tx, ty, tz;
static float scale;


int read_binvox(string filespec)
{
    ifstream *input = new ifstream(filespec.c_str(), ios::in | ios::binary);
    
    //
    // read header
    //
    string line;
    *input >> line;  // #binvox
    if (line.compare("#binvox") != 0) {
        cout << "Error: first line reads [" << line << "] instead of [#binvox]" << endl;
        delete input;
        return 0;
    }
    *input >> version;
    cout << "reading binvox version " << version << endl;
    
    depth = -1;
    int done = 0;
    while(input->good() && !done) {
        *input >> line;
        if (line.compare("data") == 0) done = 1;
        else if (line.compare("dim") == 0) {
            *input >> depth >> height >> width;
        }
        else if (line.compare("translate") == 0) {
            *input >> tx >> ty >> tz;
        }
        else if (line.compare("scale") == 0) {
            *input >> scale;
        }
        else {
            cout << "  unrecognized keyword [" << line << "], skipping" << endl;
            char c;
            do {  // skip until end of line
                c = input->get();
            } while(input->good() && (c != '\n'));
            
        }
    }
    if (!done) {
        cout << "  error reading header" << endl;
        return 0;
    }
    if (depth == -1) {
        cout << "  missing dimensions in header" << endl;
        return 0;
    }
    
    size = width * height * depth;
    voxelsVec.reserve(size);
    distanceFieldVec.reserve(size);
    
    //
    // read voxel data
    //
    byte value;
    byte count;
    int index = 0;
    int end_index = 0;
    int nr_voxels = 0;
    
    input->unsetf(ios::skipws);  // need to read every byte now (!)
    *input >> value;  // read the linefeed char
    
    while((end_index < size) && input->good())
    {
        *input >> value >> count;
        
        if (input->good())
        {
            end_index = index + count;
            if (end_index > size)
            {
                return 0;
            }
            
            for(int i=index; i < end_index; i++)
            {
                voxelsVec.push_back(value);
            }
            
            if (value)
            {
                nr_voxels += count;
            }
            
            index = end_index;
        }  // if file still ok
        
    }  // while
    
    input->close();
    
    originalVoxelsVec = voxelsVec;
    
    cout << "dim " << depth << " " << height << " " << width << endl;
    cout << "translate " << tx << " " << ty << " " << tz << endl;
    cout << "scale " << scale << endl;
    cout << "read " << nr_voxels << " voxels" << endl;
    
    
    return 1;
}

void carveVoxelObject()
{
    for (int x = 0; x < depth; x++)
    {
        for (int y = 0; y < depth; y++)
        {
            for (int z = 0; z < depth; z++)
            {
                if ((int)originalVoxelsVec[(x * height + y) * depth + z] == 1)
                {
                    bool onSurface = true;
                    
                    if ( (x > 0 && (int)originalVoxelsVec[(( x - 1) * height + y) * depth + z] == 1)
                        && ((x + 1) < depth && (int)originalVoxelsVec[(( x + 1) * height + y) * depth + z] == 1)
                        && (y > 0 && (int)originalVoxelsVec[(x * height + y - 1) * depth + z] == 1)
                        && ((y + 1) < depth && (int)originalVoxelsVec[(x * height + y + 1) * depth + z] == 1)
                        && (z > 0 && (int)originalVoxelsVec[(x * height + y) * depth + z - 1] == 1)
                        && ((z + 1) < depth && (int)originalVoxelsVec[(x * height + y) * depth + z + 1] == 1) )
                    {
                        onSurface = false;
                    }
                    
                    if (!onSurface)
                    {
                        voxelsVec[(x * height + y) * depth + z] = 0;
                    }
                }
                
                
            }
        }
    }
    
    cout << "Object has been carved." << endl;
}

// Brute force
void createDistanceField()
{
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            for (int z = 0; z < depth; z++)
            {
                double min = DBL_MAX;
                
                for (int xx = 0; xx < width; xx++)
                {
                    for (int yy = 0; yy < height; yy++)
                    {
                        for (int zz = 0; zz < depth; zz++)
                        {
                            if (voxelsVec[(xx * height + yy) * depth + zz] == 1)
                            {
                                double mag = (pow(x - xx, 2) + pow(y - yy, 2) + pow(z - zz, 2));
                                //double mag = (abs(x - xx) + abs(y - yy) + abs(z - zz));
                                
                                if (mag < min)
                                {
                                    min = mag;
                                }
                            }
                        }
                    }
                }
                
                distanceFieldVec.push_back(sqrt(min));
            }
        }
    }
    
    cout << "Distance field has been created." << endl;
}

void correctSignes()
{
    long index = 0;
    const long indexend = voxelsVec.size();
    while ( index != indexend )
    {
        if ( (int)originalVoxelsVec[index] == 1)
        {
            distanceFieldVec[index] *= -1;
        }
        
        index++;
    }
    
    cout << "Signes in distance field has been corrected." << endl;
}

void  distanceTransform3D()
{
    // feed distance field vector
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            for (int k = 0; k < depth; k++)
            {
                if (voxelsVec[(i * height + j) * depth + k] == 1)
                {
                    distanceFieldVec.push_back(0);
                    
                }
                else
                {
                    distanceFieldVec.push_back(DBL_MAX);
                    
                    // at to nepretece
                    distanceFieldVec.back() *= 0.5;
                }
            }
        }
    }
    
    // init
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            for (int k = 0; k < depth; k++)
            {
                if (distanceFieldVec[(i * height + j) * depth + k] != 0)
                {
                    distanceFieldVec[(i * height + j) * depth + k] = DBL_MAX;
                    
                    // at to nepretece
                    distanceFieldVec[(i * height + j) * depth + k] *= 0.5f;
                }
            }
        }
    }
    
    // forward
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            for (int k = 0; k < depth; k++)
            {
                double min = DBL_MAX;
                
                for (int l1 = -1; l1 <= 0; l1++)
                {
                    int g2;
                    if (l1 == -1)
                        g2 = 1;
                    else
                        g2 = 0;
                    
                    for (int l2 = -1; l2 <= g2; l2++)
                    {
                        int g3;
                        if (std::min(l1, l2) == -1)
                            g3 = 1;
                        else
                            g3 = 0;
                        
                        for (int l3 = -1; l3 <= g3; l3++)
                        {
                            int g4;
                            if (std::min({l1, l2, l3}) == -1)
                                g4 = 1;
                            else
                                g4 = 0;
                            
                            
                            // correction out of array
                            int ll1, ll2, ll3;
                            ll1 = l1;
                            ll2 = l2;
                            ll3 = l3;
                            
                            if ((i + ll1) < 0)
                            {
                                ll1 = 0;
                            }
                            else if ((i + ll1) >= width)
                            {
                                ll1 = width - 1 - i;
                            }
                            
                            if ((j + ll2) < 0)
                            {
                                ll2 = 0;
                            }
                            else if ((j + ll2) >= height)
                            {
                                ll2 = height - 1 - j;
                            }
                            
                            if ((k + ll3) < 0)
                            {
                                ll3 = 0;
                            }
                            else if ((k + ll3) >= depth)
                            {
                                ll3 = depth - 1 - k;
                            }
                            
                            min = std::min(min, distanceFieldVec[((i + ll1) * height + j + ll2) * depth + k + ll3] + (abs(l1) + abs(l2) + abs(l3)) );
                        }
                        
                        
                    }
                }
                
                distanceFieldVec[(i * height + j) * depth + k] = min;
                
            }
        }
    }
    
    
    // backward
    for (int i = width-1; i >= 0 ; i--)
    {
        for (int j = height-1; j >= 0; j--)
        {
            for (int k = depth-1; k >= 0; k--)
            {
                double min = DBL_MAX;
                
                for (int l1 = 0; l1 <= 1; l1++)
                {
                    int g2;
                    if (l1 == 1)
                        g2 = -1;
                    else
                        g2 = 0;
                    
                    for (int l2 = g2; l2 <= 1; l2++)
                    {
                        int g3;
                        if (std::min(l1, l2) == 1)
                            g3 = -1;
                        else
                            g3 = 0;
                        
                        for (int l3 = g3; l3 <= 1; l3++)
                        {
                            int g4;
                            if (std::min({l1, l2, l3}) == 1)
                                g4 = -1;
                            else
                                g4 = 0;
                            
                            
                            // correction out of array
                            int ll1, ll2, ll3;
                            ll1 = l1;
                            ll2 = l2;
                            ll3 = l3;
                            
                            if ((i + ll1) < 0)
                            {
                                ll1 = 0;
                            }
                            else if ((i + ll1) >= width)
                            {
                                ll1 = width - 1 - i;
                            }
                            
                            if ((j + ll2) < 0)
                            {
                                ll2 = 0;
                            }
                            else if ((j + ll2) >= height)
                            {
                                ll2 = height - 1 - j;
                            }
                            
                            if ((k + ll3) < 0)
                            {
                                ll3 = 0;
                            }
                            else if ((k + ll3) >= depth)
                            {
                                ll3 = depth - 1 - k;
                            }

                            min = std::min(min, distanceFieldVec[((i + ll1) * height + j + ll2) * depth + k + ll3] + (abs(l1) + abs(l2) + abs(l3)) );
                        }
                        
                        
                    }
                }
                
                distanceFieldVec[(i * height + j) * depth + k] = min;
                
            }
        }
    }
    
    cout << "Distance field has been created." << endl;
}

string createDistfieldFileName(string filePath)
{
    long pos = filePath.find("binvox");
    filePath = filePath.substr(0, pos);
    
    string suffix = "distfield";
    
    filePath.append(suffix);
    
    return filePath;
}

// struktura: dim, obj size, distance field data
void saveToFile(string filePath, float objSize)
{
    string fileName = createDistfieldFileName(filePath);
    
    ostringstream stream;
    
    // dim
    stream << width << ",";
    stream << objSize << ",";
    
    // data
    long index = 0;
    const long indexend = distanceFieldVec.size();
    while ( index != indexend )
    {
        stream << distanceFieldVec[index++] << ",";
    }
    
    ofstream myfile;
    myfile.open (fileName, ios::binary | ios::out);

    
    myfile << stream.str();//.c_str();
    
    myfile.close();
}



/*
// tests
void parse(char *newdecompressedDistanceField)
{
    vector<double> newDistanceFieldVec;
    
    char * pch;
    pch = strtok (newdecompressedDistanceField, ",");
    
    while (true)
    {
        pch = strtok (NULL, ",");
        
        if (pch != NULL)
        {
            newDistanceFieldVec.push_back(atof(pch));
        }
        else
        {
            break;
        }
    }
}

void loadFile()
{
    // load
    streampos size;
    char * memblock;
    
    ifstream file (fileName, ios::in|ios::binary|ios::ate);
    if (file.is_open())
    {
        size = file.tellg();
        memblock = new char [size];
        file.seekg (0, ios::beg);
        file.read (memblock, size);
        file.close();
        
        cout << "the entire file content is in memory" << endl;
        
        delete[] memblock;
    }
    else cout << "Unable to open file";
}
 */


// prints
void printVoxelGrid()
{
    cout << endl << "data:" << endl;
    long index = 0;
    const long indexend = voxelsVec.size();
    while ( index != indexend )
    {
        cout << (int)voxelsVec[index++] << " ";
        
        if (index % width == 0)
        {
            cout << endl;
        }
        
        if ((index % ( width * height)) == 0)
        {
            cout << endl;
        }
    }
}

void printDistanceField()
{
    
    cout << endl << "distance field:" << endl;
    long index = 0;
    const long indexend = distanceFieldVec.size();
    while ( index != indexend )
    {
        cout << distanceFieldVec[index++] << " ";
        
        if (index % width == 0)
        {
            cout << endl;
        }
        
        if ((index % ( width * height)) == 0)
        {
            cout << endl;
        }
    }
    
}


int main(int argc, const char * argv[])
{
    if (argc > 2)
    {
        clock_t begin = clock();
        
        string filePath = argv[1];
        float objSize = atof(argv[2]);
        
        
        read_binvox(filePath); //"/Users/jakubvlk/Craneballs/SandBoxTerrainEngine/TestingUnity/Assets/StreamingAssets/Models/buddha_100.binvox");
    //    printVoxelGrid();
        carveVoxelObject();
    //    printVoxelGrid();

    //    createDistanceField();
        distanceTransform3D();
        correctSignes();
    //    printDistanceField();
        
        saveToFile(filePath, objSize);
        
        // test
    //    loadFile();
        
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        
        cout << endl << "Elapsed time = " << elapsed_secs << "s." << endl;
    }
    else
    {
        cerr << "Missing parameters!" << endl;
        
        cout << "First param is file name, second param is size of obj model." << endl;
    }
    
    
    return 0;
}
