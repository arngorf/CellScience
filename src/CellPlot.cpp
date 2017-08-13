#include "CellPlot.hpp"

#include "ChanVese3D.hpp"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkCubeSource.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkImplicitVolume.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkOutlineFilter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProgrammableSource.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkReverseSense.h>
#include <vtkSampleFunction.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkTransform.h>
#include <vtkXMLPolyDataReader.h>

#include "Debug.hpp"
#include "UtilityFunctions.hpp"

CellPlot::CellPlot()
{

}

CellPlot::~CellPlot()
{

}

void CellPlot::InterativeCubeField(float *phaseField, int nz, int ny, int nx, float levelSet)
{
    // CubeActors
    std::vector<vtkSmartPointer<vtkActor>> cubeActorPtrs;

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {

                if (phaseField[gIndex(x,y,z,ny,nx,1)] < levelSet) {

                    vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
                    cubeSource->SetCenter((float) x, (float) y, (float) z);
                    cubeSource->Update();

                    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
                    mapper->SetInputConnection(cubeSource->GetOutputPort());

                    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
                    actor->SetMapper(mapper);
                    cubeActorPtrs.push_back(actor);
                }

            }
        }
    }

    // A renderer and render window
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);

    // An interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    for (int i = 0; i < cubeActorPtrs.size(); ++i) {
        renderer->AddActor(cubeActorPtrs[i]);
    }

    renderer->SetBackground(.3, .2, .1);

    // Render
    renderWindow->Render();

    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

    renderWindowInteractor->SetInteractorStyle( style );

    // Begin mouse interaction
    renderWindowInteractor->Start();
}

void CellPlot::InteractiveGuessMesh(float *phaseField, int nz, int ny, int nx)
{
    vtkSmartPointer<vtkPoints> inputPoints = vtkSmartPointer<vtkPoints>::New();

    for (int z = 1; z < nz; ++z) {
        for (int y = 1; y < ny; ++y) {
            for (int x = 1; x < nx; ++x) {


                float b = (float) phaseField[gIndex(x,y,z,ny,nx,1)];
                float a = (float) phaseField[gIndex(x-1,y,z,ny,nx,1)];


                if (a * b < 0) {
                    float interpolant = a/(a-b);
                    const float p[3] = {float(x - 1 + interpolant), (float) y, (float) z};
                    inputPoints->InsertNextPoint(p);
                }

                a = (float) phaseField[gIndex(x,y-1,z,ny,nx,1)];

                if (a * b < 0) {
                    float interpolant = a/(a-b);
                    const float p[3] = {(float) x, float(y - 1 + interpolant), (float) z};
                    inputPoints->InsertNextPoint(p);
                }

                a = (float) phaseField[gIndex(x,y,z-1,ny,nx,1)];

                if (a * b < 0) {
                    float interpolant = a/(a-b);
                    const float p[3] = {(float) x, (float) y, float(z - 1 + interpolant)};
                    inputPoints->InsertNextPoint(p);
                }
            }
        }
    }

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(inputPoints);

    // Construct the surface and create isosurface.
    vtkSmartPointer<vtkSurfaceReconstructionFilter> surf = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();

    surf->SetInputData(polydata);

    vtkSmartPointer<vtkContourFilter> cf = vtkSmartPointer<vtkContourFilter>::New();
    cf->SetInputConnection(surf->GetOutputPort());
    cf->SetValue(0, 0.0);

    // Sometimes the contouring algorithm can create a volume whose gradient
    // vector and ordering of polygon (using the right hand rule) are
    // inconsistent. vtkReverseSense cures this problem.
    vtkSmartPointer<vtkReverseSense> reverse = vtkSmartPointer<vtkReverseSense>::New();
    reverse->SetInputConnection(cf->GetOutputPort());
    reverse->ReverseCellsOn();
    reverse->ReverseNormalsOn();

    vtkSmartPointer<vtkPolyDataMapper> map = vtkSmartPointer<vtkPolyDataMapper>::New();
    map->SetInputConnection(reverse->GetOutputPort());
    map->ScalarVisibilityOff();

    vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
    surfaceActor->SetMapper(map);

    // Create the RenderWindow, Renderer and both Actors
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren);
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);

    // Add the actors to the renderer, set the background and size
    ren->AddActor(surfaceActor);
    ren->SetBackground(.2, .3, .4);

    renWin->Render();
    iren->Start();
}

void CellPlot::InteractiveImplicitSurface(float *phaseField, int nz, int ny, int nx)
{

    vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();

    imageImport->SetDataSpacing(1, 1, 1);
    imageImport->SetDataOrigin(0, 0, 0);
    imageImport->SetWholeExtent(0, nx+1, 0, ny+1, 0, nz+1);
    imageImport->SetDataExtentToWholeExtent();
    imageImport->SetDataScalarTypeToFloat();
    imageImport->SetImportVoidPointer(phaseField);
    imageImport->Update();

    // Make implicit volume function from vtk image data
    vtkSmartPointer<vtkImplicitVolume> implicitVolume = vtkSmartPointer<vtkImplicitVolume>::New();
    implicitVolume->SetVolume(imageImport->GetOutput());
    //implicitVolume->SetVolume(imageData);
    // Sample the function
    vtkSmartPointer<vtkSampleFunction> sample = vtkSmartPointer<vtkSampleFunction>::New();
    sample->SetSampleDimensions(3*(nx+1), 3*(ny+1), 3*(nz+1));

    sample->SetImplicitFunction(implicitVolume);

    sample->SetModelBounds(0, nx+1, 0, ny+1, 0, nz+1);

    // Create the 0 isosurface
    vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
    contours->SetInputConnection(sample->GetOutputPort());
    contours->GenerateValues(1, 1, 1);

    // Map the contours to graphical primitives
    vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    contourMapper->SetInputConnection(contours->GetOutputPort());
    contourMapper->ScalarVisibilityOff();

    // Create an actor for the contours
    vtkSmartPointer<vtkActor> contourActor = vtkSmartPointer<vtkActor>::New();
    contourActor->SetMapper(contourMapper);

    // Visualize
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    renderer->AddActor(contourActor);
    renderer->SetBackground(0.2, 0.3, 0.4);

    renderWindow->Render();
    interactor->Start();
}

void CellPlot::MultiObjectInteractiveImplicitSurface(std::vector<CellObject*> objects)
{
    // Visualize
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

    for (int j = 0; j < objects.size(); ++j)
    {
        CellObject *object = objects[j];

        int beginX, endX, beginY, endY, beginZ, endZ;
        object->getBounds(beginX, endX, beginY, endY, beginZ, endZ);

        int nx = endX - beginX + 3;
        int ny = endY - beginY + 3;
        int nz = endZ - beginZ + 3;

        Debug::Info("CellPlot::MultiObjectInteractiveImplicitSurface: "
                    "dimensions: " + STR(nx) + ", " + STR(ny) + ", "
                    + STR(nz));

        float *phi = object->getPhiPtr();

        if (phi == NULL)
        {
            Debug::Error("Cellplot::MultiObjectInteractiveImplicitSurface: "
                         "(float*) phi of (CellObject*) object is NULL");

            return;
        }

        vtkSmartPointer<vtkImageImport> imageImport = vtkSmartPointer<vtkImageImport>::New();

        imageImport->SetDataSpacing(1, 1, 1);
        imageImport->SetDataOrigin(0, 0, 0);
        //imageImport->SetWholeExtent(0, nx+1, 0, ny+1, 0, nz+1);
        imageImport->SetWholeExtent(0, nx-1, 0, ny-1, 0, nz-1);
        imageImport->SetDataExtentToWholeExtent();
        imageImport->SetDataScalarTypeToFloat();
        imageImport->SetImportVoidPointer(phi);
        imageImport->Update();

        // Make implicit volume function from vtk image data
        vtkSmartPointer<vtkImplicitVolume> implicitVolume = vtkSmartPointer<vtkImplicitVolume>::New();
        implicitVolume->SetVolume(imageImport->GetOutput());
        //implicitVolume->SetVolume(imageData);
        // Sample the function
        vtkSmartPointer<vtkSampleFunction> sample = vtkSmartPointer<vtkSampleFunction>::New();
        //sample->SetSampleDimensions(3*(nx+1), 3*(ny+1), 3*(nz+1));
        sample->SetSampleDimensions((nx+1), (ny+1), (nz+1));
        //sample->SetSampleDimensions(3*(nx-1), 3*(ny-1), 3*(nz-1));

        sample->SetImplicitFunction(implicitVolume);

        //sample->SetModelBounds(0, nx+1, 0, ny+1, 0, nz+1);
        sample->SetModelBounds(0, nx-1, 0, ny-1, 0, nz-1);

        // Create the 0 isosurface
        vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
        contours->SetInputConnection(sample->GetOutputPort());
        contours->GenerateValues(1, 1, 1);

        // Map the contours to graphical primitives
        vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        contourMapper->SetInputConnection(contours->GetOutputPort());
        contourMapper->ScalarVisibilityOff();

        // Create an actor for the contours
        vtkSmartPointer<vtkActor> contourActor = vtkSmartPointer<vtkActor>::New();
        contourActor->SetMapper(contourMapper);

        // Translate to proper position. The level of abstraction is... deep..
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->PostMultiply();
        transform->Translate(beginX, beginY, beginZ);
        contourActor->SetUserTransform(transform);

        // color the actor
        switch (object->getType())
        {
        case COT_UNCATEGORIZED:
            contourActor->GetProperty()->SetOpacity(0.5);
            break;
        case COT_MEMBRANE:
            contourActor->GetProperty()->SetColor(0,0.5,0.5);
            contourActor->GetProperty()->SetOpacity(0.15);
            break;
        case COT_VESICLE:
            contourActor->GetProperty()->SetColor(0.15,0.2,0.65);
            contourActor->GetProperty()->SetOpacity(0.3);
            break;
        case COT_FILAMENT:
            contourActor->GetProperty()->SetColor(1.0,0.0,1.0);
            contourActor->GetProperty()->SetOpacity(0.5);
            break;
        case COT_SYNAPSE:
            contourActor->GetProperty()->SetColor(0.7,0.2,0.1);
            contourActor->GetProperty()->SetOpacity(0.5);
            break;
        case COT_ENDOSOME:
            contourActor->GetProperty()->SetColor(0.46,0.86,0.2);
            contourActor->GetProperty()->SetOpacity(0.6);
            break;
        case COT_ER:
            contourActor->GetProperty()->SetColor(1.0,0.8,0.2);
            contourActor->GetProperty()->SetOpacity(0.5);
            break;
        }

        renderer->AddActor(contourActor);
    }

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

    //renderer->SetBackground(0.2, 0.3, 0.4);
    //renderer->SetBackground(0.0, 0.0, 0.0);
    renderer->SetBackground(1.0, 1.0, 1.0);

    // 1. Use a render window with alpha bits (as initial value is 0 (false)):
    renderWindow->SetAlphaBitPlanes(true);

    // 2. Force to not pick a framebuffer with a multisample buffer
    // (as initial value is 8):
    renderWindow->SetMultiSamples(0);

    // 3. Choose to use depth peeling (if supported) (initial value is 0 (false)):
    renderer->SetUseDepthPeeling(true);

    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    interactor->Start();
}