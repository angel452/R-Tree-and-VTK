#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkTexture.h>
#include <vtkPNGReader.h>
#include <vtkPlaneSource.h>

#include <sstream>
#include <string>
#include <fstream>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkVertexGlyphFilter.h>

#include <vtkAxesActor.h>
#include <vtkCaptionActor2D.h>
#include <vtkTransform.h>

#include "Rtree.h"

using namespace std;

// The ordering of the corner points on each face.
std::array<std::array<vtkIdType, 4>, 6> ordering = {{{{0, 3, 2, 1}},
                                                     {{4, 5, 6, 7}},
                                                     {{0, 1, 5, 4}},
                                                     {{1, 2, 6, 5}},
                                                     {{2, 3, 7, 6}},
                                                     {{3, 0, 4, 7}}}};

vtkNew<vtkActor> drawRegion( double x, double y, double z, double hx, double hy, double hz){
    std::array<std::array<double, 3>, 8> pts;
    pts[0] = {x,y,z};       //0,0,0
    pts[1] = {x+hx, y, z };  //1,0,0
    pts[2] = {x+hx, y+hy, z}; //1,1,0
    pts[3] = {x, y+hy, z};   //0,1,0
    pts[4] = {x, y, z+hz};   //0,0,1
    pts[5] = {x+hx,y,z+hz};   //1,0,1
    pts[6] = {x+hx,y+hy,z+hz}; //1,1,1
    pts[7] = {x,y+hy,z+hz};   //0,1,1

    vtkNew<vtkPolyData> cube;
    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> polys;
    vtkNew<vtkFloatArray> scalars;

    for (auto i = 0ul; i < pts.size(); ++i) {
        points->InsertPoint(i, pts[i].data());
        scalars->InsertTuple1(i, i);
    }
    for (auto&& i : ordering) {
        polys->InsertNextCell(vtkIdType(i.size()), i.data());
    }

    // We now assign the pieces to the vtkPolyData.
    cube->SetPoints(points);
    cube->SetPolys(polys);
    cube->GetPointData()->SetScalars(scalars);

    // Now we'll look at it.
    vtkNew<vtkPolyDataMapper> cubeMapper;
    cubeMapper->SetInputData(cube);
    cubeMapper->SetScalarRange(cube->GetScalarRange());
    vtkNew<vtkActor> cubeActor;
    cubeActor->SetMapper(cubeMapper);
    cubeActor->GetProperty()->SetRepresentationToWireframe();

    return cubeActor;
}

vtkNew<vtkActor> drawFigure( const string &filename, double x, double y, double z, double hx, double hy){
    vtkNew<vtkPNGReader> PNGReader;

    if (!PNGReader->CanReadFile(filename.c_str())) {
    cerr << "CanReadFile failed for " << filename.c_str() << "\n";
    exit(EXIT_FAILURE);
    }

    PNGReader->SetFileName(filename.c_str());
    PNGReader->Update();

    // objeto textura
    vtkNew<vtkTexture> texture;
    texture->SetInputConnection(PNGReader->GetOutputPort());


    // objeto plano
    vtkNew<vtkPlaneSource> planeSource;
    /*
    planeSource->SetOrigin(0, 0, 0);
    planeSource->SetPoint1(hx,0,0);
    planeSource->SetPoint2(0,hy,0);
    planeSource->SetCenter(x, y, z);
    */
    planeSource->SetOrigin(0, 0, 0);
    planeSource->SetPoint1(hx,0,0);
    planeSource->SetPoint2(0,hy,0);
    planeSource->SetCenter( x+(hx/2) , y+(hy/2) , z );

    // Esta linea lo hace un poco mas divertWido
    //planeSource->SetNormal( (rand()%100)/100.0, (rand()%100)/100.0, (rand()%100)/100.0 );
    planeSource->Update();

    // mapea y asigna textura
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(planeSource->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->SetTexture(texture);

    return actor;
}

int main(int argc, char const *argv[])
{
    //################### RTREE ##############################
    const int ndimensiones = 3; // attack - defense - speed
    int M = 50; // hijos por nodo

    cout << endl << "######################################## " << endl;
    cout << "               R - TRE " << endl;
    cout << "######################################## " << endl;
    RTree<Pokemon<ndimensiones>, ndimensiones> test1(M);
     //#########################################################

    string filename = "pokemon.csv";
    //string filename = "pokemon2.csv";

    ifstream in(filename);
    int num;
    string name;
    float attack, defense, speed, height_m, hp, weight_kg_std;
    vtkNew<vtkNamedColors> colors;
    // The usual rendering stuff.
    vtkNew<vtkCamera> camera;
    camera->SetPosition(1, 1, 1);
    camera->SetFocalPoint(0, 0, 0);

    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renWin;
    renWin->AddRenderer(renderer);
    renWin->SetWindowName("DrawPokemon");

    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetRenderWindow(renWin);


    while( in >> num >> name >> attack >> defense >> speed >> height_m >> hp >> weight_kg_std  ) {
        vtkNew<vtkActor> image = drawFigure( "pokemon_png/" + to_string(num)+".png",
                                              attack, defense, speed,
                                              height_m, hp  );
        Record<Pokemon<ndimensiones>, ndimensiones> pokemonAux(attack,defense,speed,height_m, hp, weight_kg_std);
        test1.insert(pokemonAux);
        renderer->AddActor(image);
    }

    cout << endl << "##########################################################################" << endl;
    cout << "El txt se leyo correctamente y se insertaron todos los puntos" << endl;
    cout << "##########################################################################" << endl;

    // ################### DIBUJAR REGIONES #########################
    // test
    /*
    vtkNew<vtkActor> auxActor = drawRegion(2,2,0, 2, 2, 0);
    renderer->AddActor(auxActor);

    vtkNew<vtkActor> auxActor2 = drawRegion(4,4,0, 3, 3, 0);
    renderer->AddActor(auxActor2);

    vtkNew<vtkActor> auxActor3 = drawRegion(5,2,0, 2, 2, 0);
    renderer->AddActor(auxActor3);
     */
    // fin Test
    /*
    1 Bulbasaur 2 2 0 2 2 0
    2 Ivysaur 4 4 0 3 3 0
    3 Venusaur 5 1 0 3 2 0
     */
    cout << endl << "Get all bottomLeft and upperRight to draw..." << endl;
    test1.getAllSizes();
    cout << "-----> Imprimiendo" << endl;
    for(int i = 0; i < allSizes.size(); i++){
        float bottomLeftX = allSizes[i][0];
        float bottomLeftY = allSizes[i][1];
        float bottomLeftZ = allSizes[i][2];
        float distX = allSizes[i][3];
        float distY = allSizes[i][4];
        float distZ = allSizes[i][5];

        //cout <<  endl << "##############################################" << endl;
        //cout << "Regiones resultantes: " << endl;
        //cout << bottomLeftX << " - " << bottomLeftY << " - " << bottomLeftZ << " - " << distX << " - " << distY << " - " << distZ << endl;
        vtkNew<vtkActor> auxActor = drawRegion(bottomLeftX,bottomLeftY,bottomLeftZ, distX, distY, distZ);
        renderer->AddActor(auxActor);
    }

    // ################ Para mostrar los ejes x, y, z ##################
    vtkNew<vtkTransform> transform;
    transform->Translate(1.0, 0.0, 0.0);
    vtkNew<vtkAxesActor> axes;
    // The axes are positioned with a user transform
    axes->SetUserTransform(transform);
    renderer->AddActor(axes);
    // #################################################################

    renderer->SetActiveCamera(camera);
    renderer->ResetCamera();
    renderer->SetBackground(colors->GetColor3d("Cornsilk").GetData());

    renWin->SetSize(600, 600);

    // interact with data
    renWin->Render();
    iren->Start();

    return 0;
}
