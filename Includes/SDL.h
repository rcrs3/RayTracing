#ifndef RT_SDLREADER_H
#define RT_SDLREADER_H


class SDL {
    std::string output;
    T3* eye;
    Ortho* ortho;
    Size* size;
    T3* background;
    double ambient;
    std::vector<Light*> lights;
    bool supersampling;
    double depth;
    std::vector<Object*> objects;
    void read();    //Method to read the file and fill the attribute fields

public:
    SDL(std::string path);
    ~SDL();
    std::string getOutput();
    T3* getEye();
    Ortho* getOrtho();
    Size* getSize();
    T3* getBackground();
    double getAmbient();
    std::vector<Light*> getLights();
    bool getSuperSampling();
    double getDepth();
    std::vector<Object*> getObjects();
};


#endif