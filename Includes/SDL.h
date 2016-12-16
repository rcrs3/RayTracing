using namespace std;

#ifndef RT_SDLREADER_H
#define RT_SDLREADER_H


class SDL {
    string output;
    T3 eye;
    Ortho ortho;
    Size* size;
    T3 background;
    double ambient;
    vector<Light> lights;
    bool supersampling;
    double depth;
    vector<Object> objects;
    void read();    //Method to read the file and fill the attribute fields

public:
    SDL(string path);
    ~SDL();
    string getOutput();
    T3 getEye();
    Ortho getOrtho();
    Size getSize();
    T3 getBackground();
    double getAmbient();
    vector<Light> getLights();
    bool getSuperSampling();
    double getDepth();
    vector<Object> getObjects();
};


#endif