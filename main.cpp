#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

//延时函数
#include <ctime> 
#include<iostream>
#include<cmath>
using namespace std;


void   Delay(int   time)//time*1000为秒数 
{
	clock_t   now = clock();

	while (clock() - now < time);
}


//圆锥与正方体求减

//#include<vtkSmartPointer.h>
//#include<vtkConeSource.h>
//#include<vtkImplicitPolyDataDistance.h>
//#include<vtkPointData.h>
//#include<vtkUnstructuredGrid.h>
//#include<vtkFloatArray.h>
//#include<vtkCellType.h>
//#include<vtkRectilinearGrid.h>
//#include<vtkTableBasedClipDataSet.h>
//#include<vtkPolyDataMapper.h>
//#include<vtkProperty.h>
//#include<vtkActor.h>
//#include<vtkCamera.h>
//#include<vtkRectilinearGridGeometryFilter.h>
//#include<vtkDataSetMapper.h>
//#include<vtkRenderer.h>
//#include<vtkRenderWindow.h>
//#include<vtkRenderWindowInteractor.h>
//#include<vtkCellTypes.h>
//
//#include<map>
//
//int main(int, char* [])
//{
//	auto cone =
//		vtkSmartPointer<vtkConeSource>::New();
//	cone->SetResolution(50);
//	cone->SetDirection(0, 0, -1);
//	cone->SetHeight(3.0);
//	cone->CappingOn();
//	cone->Update();
//
//	auto implicitPolyDataDistance =
//		vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
//	implicitPolyDataDistance->SetInput(cone->GetOutput());
//
//	unsigned int dimension = 51;
//	auto xCoords =
//		vtkSmartPointer<vtkFloatArray>::New();
//	for (unsigned int i = 0; i < dimension; ++i)
//	{
//		xCoords->InsertNextValue(-1 + i * 2.0 / static_cast<float>(dimension - 1));
//	}
//	auto yCoords =
//		vtkSmartPointer<vtkFloatArray>::New();
//	for (unsigned int i = 0; i < dimension; ++i)
//	{
//		yCoords->InsertNextValue(-1 + i * 2.0 / static_cast<float>(dimension - 1));
//	}
//	auto zCoords =
//		vtkSmartPointer<vtkFloatArray>::New();
//	for (unsigned int i = 0; i < dimension; ++i)
//	{
//		zCoords->InsertNextValue(-1 + i * 2.0 / static_cast<float>(dimension - 1));
//	}
//
//	auto rgrid =
//		vtkSmartPointer<vtkRectilinearGrid>::New();
//	rgrid->SetDimensions(
//		xCoords->GetNumberOfTuples(),
//		xCoords->GetNumberOfTuples(),
//		xCoords->GetNumberOfTuples());
//	rgrid->SetXCoordinates(xCoords);
//	rgrid->SetYCoordinates(yCoords);
//	rgrid->SetZCoordinates(zCoords);
//
//	auto signedDistances =
//		vtkSmartPointer<vtkFloatArray>::New();
//	signedDistances->SetNumberOfComponents(1);
//	signedDistances->SetName("SignedDistances");
//
//	for (vtkIdType pointId = 0; pointId < rgrid->GetNumberOfPoints(); ++pointId)
//	{
//		double p[3];
//		rgrid->GetPoint(pointId, p);
//		double signedDistance = implicitPolyDataDistance->EvaluateFunction(p);
//		signedDistances->InsertNextValue(signedDistance);
//	}
//
//	rgrid->GetPointData()->SetScalars(signedDistances);
//
//	auto clipper =
//		vtkSmartPointer<vtkTableBasedClipDataSet>::New();
//	clipper->SetInputData(rgrid);
//	clipper->InsideOutOn();
//	clipper->SetValue(0.0);
//	clipper->GenerateClippedOutputOn();
//	clipper->Update();
//
//	auto coneMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	coneMapper->SetInputConnection(cone->GetOutputPort());
//	auto coneActor =
//		vtkSmartPointer<vtkActor>::New();
//	coneActor->SetMapper(coneMapper);
//
//	auto geometryFilter =
//		vtkSmartPointer<vtkRectilinearGridGeometryFilter>::New();
//	geometryFilter->SetInputData(rgrid);
//	geometryFilter->SetExtent(0, dimension, 0, dimension, dimension / 2, dimension / 2);
//	geometryFilter->Update();
//
//	auto rgridMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	rgridMapper->SetInputConnection(geometryFilter->GetOutputPort());
//	rgridMapper->SetScalarRange(rgrid->GetPointData()->GetArray("SignedDistances")->GetRange());
//
//	auto wireActor =
//		vtkSmartPointer<vtkActor>::New();
//	wireActor->SetMapper(rgridMapper);
//	wireActor->GetProperty()->SetRepresentationToWireframe();
//
//	auto clipperMapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	clipperMapper->SetInputConnection(clipper->GetOutputPort());
//	clipperMapper->ScalarVisibilityOff();
//
//	auto clipperOutsideMapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	clipperOutsideMapper->SetInputConnection(clipper->GetOutputPort(1));
//	clipperOutsideMapper->ScalarVisibilityOff();
//
//	auto clipperActor =
//		vtkSmartPointer<vtkActor>::New();
//	clipperActor->SetMapper(clipperMapper);
//	clipperActor->GetProperty()->SetColor(0.89, 0.81, 0.34);
//	
//	auto clipperOutsideActor =
//		vtkSmartPointer<vtkActor>::New();
//	clipperOutsideActor->SetMapper(clipperOutsideMapper);
//	clipperOutsideActor->GetProperty()->SetColor(0.89, 0.81, 0.34);
//
//	double leftViewport[4] = { 0.0,0.0,0.5,1.0 };
//	auto leftRenderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	leftRenderer->SetViewPoint(leftViewport);
//	leftRenderer->SetBackground(.4, .5, .6);
//	leftRenderer->UseHiddenLineRemovalOn();
//
//	double rightViewport[4] = { 0.5,0.0,1.0,1.0 };
//	auto rightRenderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	rightRenderer->SetViewPoint(leftViewport);
//	rightRenderer->SetBackground(.5, .6, .6);
//	rightRenderer->UseHiddenLineRemovalOn();
//
//	leftRenderer->AddActor(wireActor);
//	leftRenderer->AddActor(clipperActor);
//	rightRenderer->AddActor(clipperOutsideActor);
//
//	auto renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(leftRenderer);
//	renderWindow->AddRenderer(rightRenderer);
//
//	auto interactor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	interactor->SetRenderWindow(renderWindow);
//
//	leftRenderer->GetActiveCamera()->SetPosition(0, -1, 0);
//	leftRenderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
//	leftRenderer->GetActiveCamera()->SetViewUp(0, 0, 1);
//	leftRenderer->GetActiveCamera()->Azimuth(30);
//	leftRenderer->GetActiveCamera()->Elevation(30);
//	leftRenderer->ResetCamera();
//	rightRenderer->SetActiveCamera(leftRenderer->GetActiveCamera());
//
//	renderWindow->Render();
//	interactor->Start();
//
//	vtkIdType numberOfCells = clipper->GetOutput()->GetNumberOfCells();
//	std::cout << "------------------------" << std::endl;
//	std::cout << "The clipped dataset(inside)contains a " << std::endl
//		<< clipper->GetOutput()->GetClassName()
//		<< "that has " << numberOfCells << "cells " << std::endl;
//	typedef std::map<int, int>CellContainer;
//	CellContainer cellMap;
//	for (vtkIdType i = 0; i < numberOfCells; ++i)
//	{
//		cellMap[clipper->GetOutput()->GetCellType(i)]++;
//	}
//
//	for (auto c : cellMap)
//	{
//		std::cout << "\tCell type "
//			<< vtkCellTypes::GetClassNameFromTypeId(c.first)
//			<< "occurs " << c.second << " times." << std::endl;
//	}
//
//	numberOfCells = clipper->GetClippedOutput()->GetNumberOfCells();
//	std::cout << "------------------------" << std::endl;
//	std::cout << "The clipped dataset(outside) contains a " << std::endl
//		<< clipper->GetClippedOutput()->GetClassName()
//		<< " that has " << numberOfCells << " cells" << std::endl;
//	typedef std::map<int, int> OutsideCellContainer;
//	CellContainer outsideCellMap;
//	for (vtkIdType i = 0; i < numberOfCells; i++)
//	{
//		outsideCellMap[clipper->GetClippedOutput()->GetCellType(i)]++;
//	}
//
//	for (auto c : outsideCellMap)
//	{
//		std::cout << "\tCell type "
//			<< vtkCellTypes::GetClassNameFromTypeId(c.first)
//			<< " occurs " << c.second << " times." << std::endl;
//	}
//
//	return EXIT_SUCCESS;
//}


/*-----------------------四窗口显示---------------------------*/
//#include <vtkConeSource.h>
//#include <vtkCubeSource.h>
//#include <vtkCylinderSource.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkActor.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//
//int main()
//{
//    vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
//    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
//    vtkSmartPointer<vtkCylinderSource> cylinder = vtkSmartPointer<vtkCylinderSource>::New();
//    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
//
//    vtkSmartPointer<vtkPolyDataMapper> coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    coneMapper->SetInputConnection(cone->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper> cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cubeMapper->SetInputConnection(cube->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper> cylinderMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    sphereMapper->SetInputConnection(sphere->GetOutputPort());
//
//    vtkSmartPointer<vtkActor> coneActor = vtkSmartPointer<vtkActor>::New();
//    coneActor->SetMapper(coneMapper);
//    vtkSmartPointer<vtkActor> cubeActor = vtkSmartPointer<vtkActor>::New();
//    cubeActor->SetMapper(cubeMapper);
//    vtkSmartPointer<vtkActor> cylinderActor = vtkSmartPointer<vtkActor>::New();
//    cylinderActor->SetMapper(cylinderMapper);
//    vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
//    sphereActor->SetMapper(sphereMapper);
//
//    vtkSmartPointer<vtkRenderer> renderer1 = vtkSmartPointer<vtkRenderer>::New();
//    renderer1->AddActor(coneActor);
//    renderer1->SetBackground(1.0, 0.3, 0.2);
//    renderer1->SetViewport(0.0, 0.0, 0.5, 0.5);
//    vtkSmartPointer<vtkRenderer> renderer2 = vtkSmartPointer<vtkRenderer>::New();
//    renderer2->AddActor(cubeActor);
//    renderer2->SetBackground(0.2, 1.0, 0.6);
//    renderer2->SetViewport(0.5, 0.0, 1.0, 0.5);
//    vtkSmartPointer<vtkRenderer> renderer3 = vtkSmartPointer<vtkRenderer>::New();
//    renderer3->AddActor(cylinderActor);
//    renderer3->SetBackground(0.2, 0.5, 1.0);
//    renderer3->SetViewport(0.0, 0.5, 0.5, 1.0);
//    vtkSmartPointer<vtkRenderer> renderer4 = vtkSmartPointer<vtkRenderer>::New();
//    renderer4->AddActor(sphereActor);
//    renderer4->SetBackground(1.0, 1.0, 0.3);
//    renderer4->SetViewport(0.5, 0.5, 1.0, 1.0);
//
//    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
//    renWin->AddRenderer(renderer1);
//    renWin->AddRenderer(renderer2);
//    renWin->AddRenderer(renderer3);
//    renWin->AddRenderer(renderer4);
//    renWin->SetSize(640, 480);
//    renWin->Render();
//    renWin->SetWindowName("ViewFour");
//
//    vtkSmartPointer<vtkRenderWindowInteractor> interactor =
//        vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    interactor->SetRenderWindow(renWin);
//
//    renWin->Render();
//    interactor->Initialize();
//    interactor->Start();
//
//    return EXIT_SUCCESS;
//}

/*-------------------圆柱和半个圆柱裁剪球体-------------------*/
//#include "vtkActor.h"
//#include "vtkRenderer.h"
//#include "vtkRenderWindow.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkCylinder.h"
//#include "vtkPlane.h"
//#include "vtkImplicitBoolean.h"
//#include "vtkPolyDataMapper.h"
//#include "vtkSphereSource.h"
//#include "vtkProperty.h"
//#include "vtkClipPolyData.h"
//#include "vtkTransformPolyDataFilter.h"
//#include "vtkTransform.h"
//#include "vtkInteractorStyleTrackballCamera.h"
//
//int main()
//{
//    vtkSphereSource* sphere = vtkSphereSource::New();
//    sphere->SetCenter(0, 0, 0);
//    sphere->SetRadius(10);
//    sphere->SetThetaResolution(40);
//    sphere->SetPhiResolution(40);
//
//    vtkCylinder* cylinder = vtkCylinder::New();//圆柱
//    cylinder->SetCenter(0, 0, 0);
//    cylinder->SetRadius(3);
//
//    vtkPlane* vPlane = vtkPlane::New();//横截面
//    vPlane->SetOrigin(0, 0, 0);
//    vPlane->SetNormal(0, -1, 0);
//
//    vtkImplicitBoolean* cuted_cylinder = vtkImplicitBoolean::New();
//    cuted_cylinder->SetOperationTypeToIntersection();
//    cuted_cylinder->AddFunction(cylinder);
//    cuted_cylinder->AddFunction(vPlane);
//
//    vtkClipPolyData* clipper = vtkClipPolyData::New();
//    clipper->SetInputConnection(sphere->GetOutputPort());
//    clipper->SetClipFunction(cylinder);
//    clipper->GenerateClipScalarsOn();
//    clipper->GenerateClippedOutputOn();
//    clipper->SetValue(0.5);
//
//    vtkTransform* transform = vtkTransform::New();
//    transform->Translate(7, 0, 0);
//    vtkTransformPolyDataFilter* filter = vtkTransformPolyDataFilter::New();
//    filter->SetInputConnection(clipper->GetOutputPort());
//    filter->SetTransform(transform);
//    vtkClipPolyData* clipper2 = vtkClipPolyData::New();
//    clipper2->SetInputConnection(filter->GetOutputPort());
//    clipper2->SetClipFunction(cuted_cylinder);
//    clipper2->GenerateClipScalarsOn();
//    clipper2->GenerateClippedOutputOn();
//    clipper2->SetValue(0.5);
//
//    vtkPolyDataMapper* map = vtkPolyDataMapper::New();
//    map->SetInputConnection(clipper2->GetOutputPort());
//    map->ScalarVisibilityOff();
//
//    vtkActor* actor = vtkActor::New();
//    actor->SetMapper(map);
//    actor->GetProperty()->SetColor(0, 1, 1);
//
//    actor->RotateX(40);
//
//    vtkRenderer* ren = vtkRenderer::New();
//    vtkRenderWindow* renWin = vtkRenderWindow::New();
//    renWin->AddRenderer(ren);
//
//    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
//    iren->SetRenderWindow(renWin);
//
//    ren->AddActor(actor);
//    ren->SetBackground(1, 1, 1);
//    renWin->SetSize(450, 450);
//
//    vtkInteractorStyleTrackballCamera* style = vtkInteractorStyleTrackballCamera::New();
//    iren->SetInteractorStyle(style);
//
//    iren->Initialize();
//    renWin->Render();
//
//    iren->Start();
//}


//导入stl文件
/*
#include "vtkCamera.h"
#include "vtkGenericRenderWindowInteractor.h"
#include "vtkInteractorStyleJoystickCamera.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLODActor.h"
#include "vtkLight.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPropPicker.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSTLReader.h"
#include "vtkShrinkPolyData.h"

int main(int argc, char* argv[]) {
	//创建绘制器对象
	vtkRenderer* ren1 = vtkRenderer::New();
	//设置相机
	ren1->GetActiveCamera()->SetClippingRange(0.294421, 29.4421);
	ren1->GetActiveCamera()->SetDistance(7.94348);
	ren1->GetActiveCamera()->SetFocalPoint(-66.9367, -49.4539, 258.453);
	ren1->GetActiveCamera()->SetPosition(-67.8091, -57.3489, 258.377);
	ren1->GetActiveCamera()->SetViewAngle(20);
	ren1->GetActiveCamera()->SetViewUp(-0.82718, 0.0860684, 0.555306);
	ren1->GetActiveCamera()->SetParallelProjection(0);
	ren1->GetActiveCamera()->SetUseHorizontalViewAngle(0);
	ren1->SetBackground(0.1, 0.2, 0.4);
	ren1->SetLightFollowCamera(1);
	//创建绘制窗口
	vtkRenderWindow* renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren1);
	renWin->SetSize(1134, 624);
	//创建交互器
	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	iren->SetLightFollowCamera(1);
	//读源对象读取stl数据文件
	vtkSTLReader* part = vtkSTLReader::New();
	part->SetOutput(part->GetOutput());
	part->SetFileName("zhuangpei.stl");
	//创建过滤器对象，该对象将输入数据集的每个单元向单元质心收缩
	//将会导致相邻单元之间出现裂缝
	vtkShrinkPolyData* shrink = vtkShrinkPolyData::New();
	//将源对象和过滤器连接
	shrink->SetInputData((vtkPolyData*)part->GetOutput());
	//设置收缩系数，如果为1，不收缩
	shrink->SetShrinkFactor(0.9);
	//创建映射器对象
	vtkPolyDataMapper* partMapper = vtkPolyDataMapper::New();
	partMapper->SetInputData((vtkPolyData*)shrink->GetOutput());
	partMapper->SetNumberOfPieces(1);
	partMapper->SetScalarRange(0, 1);
	partMapper->SetColorMode(0);
	partMapper->SetResolveCoincidentTopology(0);
	partMapper->SetScalarMode(0);
	//	partMapper->SetImmediateModeRendering(0);
	partMapper->SetScalarVisibility(1);
	partMapper->SetUseLookupTableScalarRange(0);
	//创建Props对象(Actor)
	vtkLODActor* partActor = vtkLODActor::New();
	partActor->SetMapper(partMapper);
	partActor->GetProperty()->SetAmbientColor(0.8275, 0.8275, 0.8275);
	partActor->GetProperty()->SetColor(0.8275, 0.8275, 0.8275);
	partActor->GetProperty()->SetDiffuseColor(0.8275, 0.8275, 0.8275);
	partActor->GetProperty()->SetOpacity(1);
	partActor->GetProperty()->SetInterpolation(1);
	partActor->GetProperty()->SetRepresentation(2);
	partActor->GetProperty()->SetBackfaceCulling(0);
	partActor->GetProperty()->SetEdgeVisibility(0);
	partActor->GetProperty()->SetFrontfaceCulling(0);
	partActor->SetOrigin(0, 0, 0);
	partActor->SetPosition(0, 0, 0);
	partActor->SetScale(1, 1, 1);
	partActor->SetVisibility(1);
	//将Actor对象添加到绘制器中
	ren1->AddActor(partActor);
	//绘制
	ren1->ResetCamera();
	ren1->ResetCameraClippingRange();
	renWin->Render();
	iren->Initialize();
	iren->Start();
	//删除对象
	iren->Delete();
	part->Delete();
	partActor->Delete();
	partMapper->Delete();
	ren1->Delete();
	renWin->Delete();
	shrink->Delete();
	return 0;
}*/

////读取stl文件

//#include <vtkSmartPointer.h>
//#include <vtkRendererCollection.h>
//#include <vtkPointPicker.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkActor.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkObjectFactory.h>
//#include <vtkProperty.h>
//#include <vtkSTLReader.h>
//#include <vtkCylinderSource.h>
//#include <vtkCubeAxesActor.h>
//#include <vtkAxesActor.h>
//#include <vtkOrientationMarkerWidget.h>
//#include <vtkAutoInit.h> 
//
//int main(int, char* [])
//{
//	 //Read a stl file.
//	vtkSmartPointer<vtkPolyData> input1 = vtkSmartPointer<vtkPolyData>::New();
//	vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
//	reader1->SetFileName("zhuangpei.stl");
//	reader1->Update();
//	input1->DeepCopy(reader1->GetOutput());
//
//	 //Create a mapper and actor
//	vtkSmartPointer<vtkPolyDataMapper> mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(reader1->GetOutputPort());
//	vtkSmartPointer<vtkActor> actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetOpacity(1);
//
//	 //Create a renderer, render window, and interactor
//	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->Render();
//	renderWindow->SetWindowName("Test");
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();
//
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetPicker(pointPicker);
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	 //Add the actor to the scene
//	renderer->AddActor(actor);
//	renderer->SetBackground(0, 0, 0);
//
//	 //Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//流线显示，找不到vtkStreamLine

//#include <vtkVersion.h>
//#include <vtkSmartPointer.h>
//#include <vtkActor.h>
//#if VTK_MAJOR_VERSION <= 5
//#include <vtkPLOT3DReader.h>
//#else
//#include <vtkMultiBlockPLOT3DReader.h>
//#include <vtkMultiBlockDataSet.h>
//#endif
//#include <vtkPlaneSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkStreamLine.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkStructuredGridOutlineFilter.h>
//#include <vtkProperty.h>
//
//int main(int argc, char* argv[])
//{
//	if (argc < 3)
//	{
//		std::cerr << "Required arguments: xyzFile qFile" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	std::string xyzFile = argv[1]; // "combxyz.bin";
//	std::string qFile = argv[2]; // "combq.bin";
//
//#if VTK_MAJOR_VERSION <= 5
//	vtkSmartPointer<vtkPLOT3DReader> pl3d =
//		vtkSmartPointer<vtkPLOT3DReader>::New();
//#else
//	vtkSmartPointer<vtkMultiBlockPLOT3DReader> pl3d =
//		vtkSmartPointer<vtkMultiBlockPLOT3DReader>::New();
//#endif
//	pl3d->SetXYZFileName(xyzFile.c_str());
//	pl3d->SetQFileName(qFile.c_str());
//	pl3d->SetScalarFunctionNumber(100);
//	pl3d->SetVectorFunctionNumber(202);
//	pl3d->Update();
//
//	// Source of the streamlines
//	vtkSmartPointer<vtkPlaneSource> seeds =
//		vtkSmartPointer<vtkPlaneSource>::New();
//	seeds->SetXResolution(4);
//	seeds->SetYResolution(4);
//	seeds->SetOrigin(2, -2, 26);
//	seeds->SetPoint1(2, 2, 26);
//	seeds->SetPoint2(2, -2, 32);
//
//	// Streamline itself
//	vtkSmartPointer<vtkStreamLine> streamLine =
//		vtkSmartPointer<vtkStreamLine>::New();
//#if VTK_MAJOR_VERSION <= 5
//	streamLine->SetInputConnection(pl3d->GetOutputPort());
//	streamLine->SetSource(seeds->GetOutput());
//#else
//	pl3d->Update();
//	streamLine->SetInputData(pl3d->GetOutput()->GetBlock(0));
//	streamLine->SetSourceConnection(seeds->GetOutputPort());
//#endif
//	//streamLine->SetStartPosition(2,-2,30);
//	// as alternative to the SetSource(), which can handle multiple
//	// streamlines, you can set a SINGLE streamline from
//	// SetStartPosition()
//	streamLine->SetMaximumPropagationTime(200);
//	streamLine->SetIntegrationStepLength(.2);
//	streamLine->SetStepLength(.001);
//	streamLine->SetNumberOfThreads(1);
//	streamLine->SetIntegrationDirectionToForward();
//	streamLine->VorticityOn();
//
//	vtkSmartPointer<vtkPolyDataMapper> streamLineMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	streamLineMapper->SetInputConnection(streamLine->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> streamLineActor =
//		vtkSmartPointer<vtkActor>::New();
//	streamLineActor->SetMapper(streamLineMapper);
//	streamLineActor->VisibilityOn();
//
//	// Outline-Filter for the grid
//	vtkSmartPointer<vtkStructuredGridOutlineFilter> outline =
//		vtkSmartPointer<vtkStructuredGridOutlineFilter>::New();
//#if VTK_MAJOR_VERSION <= 5
//	outline->SetInputConnection(pl3d->GetOutputPort());
//#else
//	outline->SetInputData(pl3d->GetOutput()->GetBlock(0));
//#endif
//	vtkSmartPointer<vtkPolyDataMapper> outlineMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	outlineMapper->SetInputConnection(outline->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> outlineActor =
//		vtkSmartPointer<vtkActor>::New();
//	outlineActor->SetMapper(outlineMapper);
//	outlineActor->GetProperty()->SetColor(1, 1, 1);
//
//	// Create the RenderWindow, Renderer and Actors
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	interactor->SetRenderWindow(renderWindow);
//
//	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
//		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
//	interactor->SetInteractorStyle(style);
//
//	renderer->AddActor(streamLineActor);
//	renderer->AddActor(outlineActor);
//
//	// Add the actors to the renderer, set the background and size
//	renderer->SetBackground(0.1, 0.2, 0.4);
//	renderWindow->SetSize(300, 300);
//	interactor->Initialize();
//	renderWindow->Render();
//
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}


// Demonstrate moving pieces and "snapping"

//#include <vtkSelectEnclosedPoints.h>
//#include <vtkRendererCollection.h>
//#include <vtkPointData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkTransform.h>
//#include <vtkLinearTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkObjectFactory.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkCubeSource.h>
//#include <vtkInteractorStyleTrackballActor.h>
//#include <vtkInteractorStyleSwitch.h>
//
//// Define interaction style
//class MouseInteractorStyle6 : public vtkInteractorStyleTrackballActor
//{
//public:
//	static MouseInteractorStyle6* New();
//	vtkTypeMacro(MouseInteractorStyle6, vtkInteractorStyleTrackballActor);
//
//	virtual void OnLeftButtonDown() override
//	{
//		std::cout << "Pressed left mouse button." << std::endl;
//		// Forward events
//		vtkInteractorStyleTrackballActor::OnLeftButtonDown();
//	}
//
//	virtual void OnMiddleButtonUp() override
//	{
//		//std::cout << "Pressed middle mouse button." << std::endl;
//
//		int x = this->Interactor->GetEventPosition()[0];
//		int y = this->Interactor->GetEventPosition()[1];
//		this->FindPokedRenderer(x, y);
//		this->FindPickedActor(x, y);
//
//		if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
//		{
//			std::cout << "Nothing selected." << std::endl;
//			return;
//		}
//
//		vtkSmartPointer<vtkPropCollection> actors =
//			vtkSmartPointer<vtkPropCollection>::New();
//
//		this->InteractionProp->GetActors(actors);
//		actors->InitTraversal();
//		vtkActor* actor = dynamic_cast<vtkActor*>(actors->GetNextProp());
//
//		vtkPolyData* polydata = dynamic_cast<vtkPolyData*>(actor->GetMapper()->GetInputAsDataSet());
//
//		vtkSmartPointer<vtkTransform> transform =
//			vtkSmartPointer<vtkTransform>::New();
//		transform->SetMatrix(actor->GetMatrix());
//
//		vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyData =
//			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//		transformPolyData->SetInputData(polydata);
//		transformPolyData->SetTransform(transform);
//		transformPolyData->Update();
//
//		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =
//			vtkSmartPointer<vtkSelectEnclosedPoints>::New();
//		selectEnclosedPoints->SetInputConnection(transformPolyData->GetOutputPort());
//		selectEnclosedPoints->SetSurfaceData(this->Sphere);
//		selectEnclosedPoints->Update();
//
//		vtkDataArray* insideArray = dynamic_cast<vtkDataArray*>(selectEnclosedPoints->GetOutput()->GetPointData()->GetArray("SelectedPoints"));
//
//		bool inside = false;
//		for (vtkIdType i = 0; i < insideArray->GetNumberOfTuples(); i++)
//		{
//			if (insideArray->GetComponent(i, 0) == 1)
//			{
//				inside = true;
//				break;
//			}
//		}
//
//		if (inside)
//		{
//			std::cout << "A point of the cube is inside the sphere!" << std::endl;
//			// Reset the cube to its original position
//			//this->CubeActor->GetMatrix()->Identity();
//			//this->CubeActor->SetOrigin(0,0,0);
//			this->CubeActor->SetPosition(0, 0, 0);
//			this->CubeActor->SetOrientation(0, 0, 0);
//
//			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->Render();
//			this->Interactor->GetRenderWindow()->Render();
//		}
//
//		// Release interaction
//		this->StopState();
//
//	}
//
//	virtual void OnRightButtonDown() override
//	{
//		std::cout << "Pressed right mouse button." << std::endl;
//		// Forward events
//		vtkInteractorStyleTrackballActor::OnRightButtonDown();
//	}
//
//	vtkPolyData* Sphere;
//	vtkActor* CubeActor;
//};
//vtkStandardNewMacro(MouseInteractorStyle6);
//
//int main(int, char* [])
//{
//	// Sphere
//	vtkSmartPointer<vtkSphereSource> sphereSource =
//		vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource->SetRadius(2);
//	sphereSource->Update();
//
//	vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> sphereActor =
//		vtkSmartPointer<vtkActor>::New();
//	sphereActor->SetMapper(sphereMapper);
//
//	// Cube
//	vtkSmartPointer<vtkCubeSource> cubeSource =
//		vtkSmartPointer<vtkCubeSource>::New();
//	cubeSource->SetCenter(5.0, 0.0, 0.0);
//	cubeSource->Update();
//
//	vtkSmartPointer<vtkPolyDataMapper> cubeMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> cubeActor =
//		vtkSmartPointer<vtkActor>::New();
//	cubeActor->SetMapper(cubeMapper);
//
//	// Visualize
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(sphereActor);
//	renderer->AddActor(cubeActor);
//
//	renderWindow->Render();
//
//	vtkSmartPointer<MouseInteractorStyle6> style =
//		vtkSmartPointer<MouseInteractorStyle6>::New();
//	style->Sphere = sphereSource->GetOutput();
//	style->CubeActor = cubeActor;
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


/*------------------stl上色--------------------------*/
//#include <vtkActor.h>
//#include <vtkFloatArray.h>
//#include <vtkLookupTable.h>
//#include <vtkPointData.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//#include <vtkSTLReader.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include<vtkUnstructuredGridReader.h>
//
//int main(int, char* [])
//{
//	// 加载一个J20的STL模型
//	vtkSmartPointer<vtkSTLReader> source = vtkSmartPointer<vtkSTLReader>::New();
//	source->SetFileName("zhuangpei.stl");
//	source->Update();
//
//	int numPts = source->GetOutput()->GetPoints()->GetNumberOfPoints();					// 获取模型的顶点数量
//	vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();		// 创建存储顶点属性的float数组
//	scalars->SetNumberOfValues(numPts);
//	for (int i = 0; i < numPts; ++i)		// 为属性数组中的每个元素设置标量值（这个标量值可以当作颜色值）
//		scalars->SetValue(i, i); //这个地方点和标量产生联系(第一个参数是第几个点；第二个参数就是设定的标量值)现在就假设第i个点颜色就是i
//
//		vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
//	poly->DeepCopy(source->GetOutput());
//	poly->GetPointData()->SetScalars(scalars);
//
//	// 创建颜色查找表
//	vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
//	hueLut->SetNumberOfColors(numPts);		// 指定颜色查找表中有多少种颜色
//	hueLut->SetHueRange(0.6667, 0.0);	//蓝到红渐变
//	hueLut->Build();
//
//	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputData(poly);
//	mapper->SetScalarRange(0, numPts);			// 设置标量值的范围（0-100）
//	mapper->ScalarVisibilityOn();
//	//mapper->SetColorModeToMapScalars();		// 无论变量数据是何种类型，该方法都通过查询表对标量数据进行映射
//	mapper->SetColorModeToDefault();			// 默认的映射器行为，即把unsigned char类型的标量属性数据当作颜色值，不执行隐式。对于其他类型的标量数据，将通过查询表映射。
//	mapper->SetLookupTable(hueLut);
//
//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//	renderer->GradientBackgroundOn();
//	renderer->SetBackground(1, 1, 1);
//	renderer->SetBackground2(1, 1, 1);
//
//	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderer->AddActor(actor);
//	renderWindow->SetSize(600, 600);
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return 0;
//}




////stl文件染色

//#include <vtkVersion.h>
//#include <vtkPlaneSource.h>
//#include <vtkPolyData.h>
//#include <vtkSTLReader.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkConeSource.h>
//#include <vtkCellArray.h>
//#include <vtkFloatArray.h>
//#include "vtkPointData.h"
//#include "vtkPoints.h"
//#include "vtkPolyData.h"
//#include "vtkPolyDataMapper.h"
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkLookupTable.h>
//
//using namespace std;
//
//int main(int, char* [])
//{
//
//	std::string inputFilename = "zhuangpei.stl";//stl文件路径
//
//	//读取stl文件
//	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//	reader->SetFileName(inputFilename.c_str());
//	reader->Update();
//
//	//创建mapper
//	vtkSmartPointer<vtkPolyDataMapper> mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(reader->GetOutputPort());
//
//
//	//创建对象
//	vtkPolyData* cube = vtkPolyData::New();//立方体
//	vtkPoints* points = vtkPoints::New();//顶点坐标
//	vtkCellArray* polys = vtkCellArray::New();//单元数组
//
//	//存储标量值
//	vtkFloatArray* scalars = vtkFloatArray::New();
//
//	vtkPolyData* polydata = reader->GetOutput();
//
//	//用于动态创建二维数组,用于存储顶点的索引
//	vtkIdType(*y)[3] = new vtkIdType[polydata->GetNumberOfCells()][3];
//
//	//读取细胞单元，并将坐标的索引赋值给二维数组y
//	for (int i = 0; i < polydata->GetNumberOfCells(); i++)
//	{
//		//l1=polydata->GetCell(i)->GetPointIds();
//		y[i][0] = polydata->GetCell(i)->GetPointIds()->GetId(0);
//
//		y[i][1] = polydata->GetCell(i)->GetPointIds()->GetId(1);
//
//		y[i][2] = polydata->GetCell(i)->GetPointIds()->GetId(2);
//
//	}
//
//	//存储顶点
//	for (int i = 0; i < polydata->GetNumberOfPoints(); i++)
//	{
//		double x[] = { 0,0,0 };
//		polydata->GetPoint(i, x);//获取顶点坐标
//		points->InsertPoint(i, x);//将顶点坐标插入到vtkPoints定义的points
//	}
//
//	//设定单元
//	for (int i = 0; i < polydata->GetNumberOfCells(); i++)
//	{
//		polys->InsertNextCell(3, y[i]);
//	}
//	//存储每个顶点的标量值,也就是颜色的索引值，这里暂时是以顶点的先后顺序来设定的
//	for (int i = 0; i < polydata->GetNumberOfPoints(); i++)
//	{
//		scalars->InsertTuple1(i, i);
//	}
//
//	//创建多边形数据
//	cube->SetPoints(points);
//	//设定单元类型为多边形
//	cube->SetPolys(polys);
//	//设定每个顶点的标量值
//	cube->GetPointData()->SetScalars(scalars);
//
//	points->Delete();
//	polys->Delete();
//	scalars->Delete();
//
//	//定义颜色映射表
//	vtkLookupTable* pColorTable = vtkLookupTable::New();
//
//	//设置颜色表中的颜色，下列两种方式都可以
//	/*
//	pColorTable->SetNumberOfColors(4);
//	pColorTable->SetTableValue(0,1.0,0.0,0.0,1.0);
//	pColorTable->SetTableValue(0,1.0,0.0,0.0,1.0);
//	pColorTable->SetTableValue(1,0.0,1.0,0.0,1.0);
//	pColorTable->SetTableValue(2,1.0,1.0,0.0,1.0);
//	pColorTable->SetTableValue(3,0.0,0.0,1.0,1.0);
//	*/
//	//设置颜色表中的颜色
//	pColorTable->SetNumberOfColors(256);
//	pColorTable->SetHueRange(0.67, 0.0);        //色调范围从红色到蓝色
//
//	pColorTable->Build();
//
//	//数据映射
//	vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
//	cubeMapper->SetInputData(cube);
//	cubeMapper->SetScalarRange(0, polydata->GetNumberOfPoints() - 1);
//	cubeMapper->SetLookupTable(pColorTable);
//	vtkActor* cubeActor = vtkActor::New();
//	cubeActor->SetMapper(cubeMapper);
//
//
//	// Create a renderer(渲染器), render window and interactor(渲染窗口)
//	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetSize(1300, 1300);//设置窗口大小
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the scene
//	renderer->AddActor(cubeActor);
//	renderer->SetBackground(.1, .2, .3); // Background color dark blue
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//
//
//	return EXIT_SUCCESS;
//}



//#include <vtkVersion.h>
//#include <vtkPlaneSource.h>
//#include <vtkPolyData.h>
//#include <vtkSTLReader.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkConeSource.h>
//#include <vtkCellArray.h>
//#include <vtkFloatArray.h>
//#include "vtkPointData.h"
//#include "vtkPoints.h"
//#include "vtkPolyData.h"
//#include "vtkPolyDataMapper.h"
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkLookupTable.h>
//#include<vtkPolyDataReader.h>
//#include<vtkUnstructuredGridReader.h>
//#include<vtkDataSetReader.h>
//#include<vtkStructuredPointsReader.h>
//
//using namespace std;
//
//int main(int, char* [])
//{
//
//	std::string filename = "carotid.vtk";
//
//	vtkSmartPointer<vtkRenderer > aRenderer =
//		vtkSmartPointer<vtkRenderer >::New();
//	vtkSmartPointer<vtkRenderWindow > renWin =
//		vtkSmartPointer<vtkRenderWindow >::New();
//	renWin->AddRenderer(aRenderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor > iren =
//		vtkSmartPointer<vtkRenderWindowInteractor >::New();
//	iren->SetRenderWindow(renWin);
//
//	vtkSmartPointer<vtkStructuredPointsReader > vtkReader = vtkSmartPointer<vtkStructuredPointsReader >::New();
//	vtkReader->SetFileName(filename.c_str());
//	vtkReader->Update();
//	vtkReader->GetReadAllScalars();
//
//	vtkSmartPointer<vtkPolyDataMapper > skinMapper = vtkSmartPointer<vtkPolyDataMapper >::New();
//	skinMapper->SetInputConnection(vtkReader->GetOutputPort());
//	skinMapper->ScalarVisibilityOff();
//
//	vtkSmartPointer<vtkActor > skin =
//		vtkSmartPointer<vtkActor >::New();
//	skin->SetMapper(skinMapper);
//
//	vtkSmartPointer<vtkCamera > aCamera =
//		vtkSmartPointer<vtkCamera >::New();
//	aCamera->SetViewUp(0, 0, -1);
//	aCamera->SetPosition(0, 1, 0);
//	aCamera->SetFocalPoint(0, 0, 0);
//	aCamera->ComputeViewPlaneNormal();
//	aCamera->Azimuth(30.0);
//	aCamera->Elevation(30.0);
//	aCamera->Dolly(1.5);
//
//	aRenderer->AddActor(skin);
//	aRenderer->SetActiveCamera(aCamera);
//	aRenderer->ResetCamera();
//	aRenderer->SetBackground(.2, .3, .4);
//	aRenderer->ResetCameraClippingRange();
//
//	renWin->Render();
//	iren->Initialize();
//	iren->Start();
//}


//手写文件上色

//#include "vtkActor.h"
//#include "vtkCamera.h"
//#include "vtkConeSource.h"
//#include "vtkCellArray.h"
//#include "vtkFloatArray.h"
//#include "vtkPointData.h"
//#include "vtkPoints.h"
//#include "vtkPolyData.h"
//#include "vtkPolyDataMapper.h"
//#include "vtkRenderWindow.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkRenderer.h"
//#include <vtkLookupTable.h>
//int main()
//{
//	int i;
//	//梯形的顶点坐标
//	static float x[8][3] = { {0,0,0}, {4,0,0}, {4,4,0}, {0,4,0},
//	{1,1,1}, {3,1,1}, {3,3,1}, {1,3,1} };
//	//4个顶点构成一个单元，一共6个单元
//	static vtkIdType y[6][4] = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4},
//	{1,2,6,5}, {2,3,7,6}, {3,0,4,7} };
//
//	vtkPoints* points = vtkPoints::New();
//	for (i = 0; i < 8; i++)
//		points->InsertPoint(i, x[i]);
//
//	vtkCellArray* polys = vtkCellArray::New();
//	for (i = 0; i < 6; i++)
//		polys->InsertNextCell(4, y[i]);
//	//存储标量值
//	vtkFloatArray* scalars = vtkFloatArray::New();
//	for (i = 0; i < 8; i++)
//		scalars->InsertTuple1(i, i);
//	//构建多边形数据
//	vtkPolyData* cube = vtkPolyData::New();
//	cube->SetPoints(points);
//	//设定单元的组成方式
//	cube->SetPolys(polys);
//	cube->GetPointData()->SetScalars(scalars);
//
//	//定义颜色映射表
//	vtkLookupTable* pColorTable = vtkLookupTable::New();
//	pColorTable->SetNumberOfColors(6);
//	pColorTable->SetTableValue(0, 1.0, 0.0, 1.0, 1.0);
//	pColorTable->SetTableValue(1, 0.0, 1.0, 1.0, 1.0);
//	pColorTable->SetTableValue(2, 1.0, 1.0, 1.0, 1.0);
//	pColorTable->SetTableValue(3, 1.0, 0.0, 1.0, 1.0);
//	pColorTable->SetTableValue(4, 0.0, 0.0, 1.0, 1.0);
//	pColorTable->SetTableValue(5, 1.0, 1.0, 0.0, 1.0);
//	pColorTable->Build();
//
//	//数据映射
//	vtkPolyDataMapper* cubeMapper = vtkPolyDataMapper::New();
//	//cubeMapper->SetInput(cube);
//	cubeMapper->SetInputData(cube);
//	cubeMapper->SetScalarRange(0, 7);
//	cubeMapper->SetLookupTable(pColorTable);
//	vtkActor* cubeActor = vtkActor::New();
//	cubeActor->SetMapper(cubeMapper);
//
//	vtkCamera* camera = vtkCamera::New();
//	camera->SetPosition(1, 1, 1);
//	camera->SetFocalPoint(0, 0, 0);
//
//	vtkRenderer* renderer = vtkRenderer::New();
//	vtkRenderWindow* renWin = vtkRenderWindow::New();
//	renWin->AddRenderer(renderer);
//
//	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
//	iren->SetRenderWindow(renWin);
//	renderer->AddActor(cubeActor);
//	renderer->SetActiveCamera(camera);
//	renderer->ResetCamera();
//	renderer->SetBackground(1, 1, 1);
//	renWin->SetSize(400, 400);
//	renWin->Render();
//	iren->Start();
//	//删除
//	points->Delete();
//	polys->Delete();
//	scalars->Delete();
//	cube->Delete();
//	cubeMapper->Delete();
//	cubeActor->Delete();
//	camera->Delete();
//	renderer->Delete();
//	renWin->Delete();
//	iren->Delete();
//	pColorTable->Delete();
//	return 0;
//}


//切割显示vtk文件

//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkCellTypes.h>
//#include <vtkDataSetMapper.h>
//#include <vtkLookupTable.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPlane.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkTableBasedClipDataSet.h>
//#include <vtkTransform.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//
//int main(int argc, char* argv[])
//{
//	//if (argc < 2)
//	//{
//	//	std::cout << "Usage: " << argv[0] << " filename.vtk e.g. treemesh.vtk"
//	//		<< std::endl;
//	//	return EXIT_FAILURE;
//	//}
//	//// Create the reader for the data.
//	//std::string filename = argv[1];
//	//std::cout << "Loading " << filename.c_str() << std::endl;
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//reader->SetFileName(filename.c_str());
//	reader->SetFileName("Test1.vtk");
//	reader->Update();
//
//	double bounds[6];
//	reader->GetOutput()->GetBounds(bounds);
//	double center[3];
//	reader->GetOutput()->GetCenter(center);
//
//	vtkNew<vtkNamedColors> colors;
//	vtkNew<vtkRenderer> renderer;
//	renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetSize(640, 480);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	double xnorm[3] = { -1.0, -1.0, 1.0 };
//
//	vtkNew<vtkPlane> clipPlane;
//	clipPlane->SetOrigin(reader->GetOutput()->GetCenter());
//	clipPlane->SetNormal(xnorm);
//
//	vtkNew<vtkTableBasedClipDataSet> clipper;
//	clipper->SetClipFunction(clipPlane);
//	clipper->SetInputData(reader->GetOutput());
//	clipper->SetValue(0.0);
//	clipper->GenerateClippedOutputOn();
//	clipper->Update();
//
//	vtkNew<vtkDataSetMapper> insideMapper;
//	insideMapper->SetInputData(clipper->GetOutput());
//	insideMapper->ScalarVisibilityOff();
//
//	vtkNew<vtkActor> insideActor;
//	insideActor->SetMapper(insideMapper);
//	insideActor->GetProperty()->SetDiffuseColor(
//		colors->GetColor3d("Banana").GetData());
//	insideActor->GetProperty()->SetAmbient(.3);
//	insideActor->GetProperty()->EdgeVisibilityOn();
//
//	vtkNew<vtkDataSetMapper> clippedMapper;
//	clippedMapper->SetInputData(clipper->GetClippedOutput());
//	clippedMapper->ScalarVisibilityOff();
//
//	vtkNew<vtkActor> clippedActor;
//	clippedActor->SetMapper(clippedMapper);
//	clippedActor->GetProperty()->SetDiffuseColor(
//		colors->GetColor3d("Tomato").GetData());
//	insideActor->GetProperty()->SetAmbient(.3);
//	clippedActor->GetProperty()->EdgeVisibilityOn();
//
//	// Create transforms to make a better visualization
//	vtkNew<vtkTransform> insideTransform;
//	insideTransform->Translate(-(bounds[1] - bounds[0]) * .75, 0, 0);
//	insideTransform->Translate(center[0], center[1], center[2]);
//	insideTransform->RotateY(-120.0);
//	insideTransform->Translate(-center[0], -center[1], -center[2]);
//	insideActor->SetUserTransform(insideTransform);
//
//	vtkNew<vtkTransform> clippedTransform;
//	clippedTransform->Translate((bounds[1] - bounds[0]) * .75, 0, 0);
//	clippedTransform->Translate(center[0], center[1], center[2]);
//	clippedTransform->RotateY(60.0);
//	clippedTransform->Translate(-center[0], -center[1], -center[2]);
//	clippedActor->SetUserTransform(clippedTransform);
//
//	renderer->AddViewProp(clippedActor);
//	renderer->AddViewProp(insideActor);
//
//	renderer->ResetCamera();
//	renderer->GetActiveCamera()->Dolly(1.4);
//	renderer->ResetCameraClippingRange();
//	renderWindow->Render();
//	renderWindow->SetWindowName("ClipUnstructuredGridWithPlane");
//	renderWindow->Render();
//
//	interactor->Start();
//
//	// Generate a report
//	vtkIdType numberOfCells = clipper->GetOutput()->GetNumberOfCells();
//	std::cout << "------------------------" << std::endl;
//	std::cout << "The inside dataset contains a " << std::endl
//		<< clipper->GetOutput()->GetClassName() << " that has "
//		<< numberOfCells << " cells" << std::endl;
//	typedef std::map<int, int> CellContainer;
//	CellContainer cellMap;
//	for (vtkIdType i = 0; i < numberOfCells; i++)
//	{
//		cellMap[clipper->GetOutput()->GetCellType(i)]++;
//	}
//
//	for (auto c : cellMap)
//	{
//		std::cout << "\tCell type " << vtkCellTypes::GetClassNameFromTypeId(c.first)
//			<< " occurs " << c.second << " times." << std::endl;
//	}
//
//	numberOfCells = clipper->GetClippedOutput()->GetNumberOfCells();
//	std::cout << "------------------------" << std::endl;
//	std::cout << "The clipped dataset contains a " << std::endl
//		<< clipper->GetClippedOutput()->GetClassName() << " that has "
//		<< numberOfCells << " cells" << std::endl;
//	typedef std::map<int, int> OutsideCellContainer;
//	CellContainer outsideCellMap;
//	for (vtkIdType i = 0; i < numberOfCells; i++)
//	{
//		outsideCellMap[clipper->GetClippedOutput()->GetCellType(i)]++;
//	}
//
//	for (auto c : outsideCellMap)
//	{
//		std::cout << "\tCell type " << vtkCellTypes::GetClassNameFromTypeId(c.first)
//			<< " occurs " << c.second << " times." << std::endl;
//	}
//	return EXIT_SUCCESS;
//}



//读取vtk文件

//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//
//#include <algorithm>
//#include <array>
//#include <string>
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//reader->SetFileName("mesh-paraview.vtk");
//	reader->SetFileName("1.vtk");
//	reader->Update();
//
//	//std::cout << "Loading: " << argv[1] << std::endl;
//	//auto unstructuredGrid = ReadUnstructuredGrid(std::string(argv[1]));
//	auto unstructuredGrid = reader->GetOutput();
//
//	// Visualize
//	vtkNew<vtkDataSetMapper> mapper;
//	mapper->SetInputData(unstructuredGrid);
//	mapper->ScalarVisibilityOff();
//
//	vtkNew<vtkProperty> backProp;
//	backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
//	backProp->SetSpecular(.6);
//	backProp->SetSpecularPower(30);
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->SetBackfaceProperty(backProp);
//	actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
//	actor->GetProperty()->SetSpecular(.3);
//	actor->GetProperty()->SetSpecularPower(30);
//	actor->GetProperty()->EdgeVisibilityOn();
//	actor->GetProperty()->SetRepresentationToWireframe();//显示外网格
//
//	renderer->AddActor(actor);
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}


//六视图

//#include <vtkSmartPointer.h>
//#include <vtkCameraActor.h>
//#include <vtkNamedColors.h>
//#include <vtkSphereSource.h>
//#include <vtkConeSource.h>
//#include <vtkCubeSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkCamera.h>
//#include <vtkMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//
//void ViewDirection(vtkRenderer* renderer,
//	double lookX, double lookY, double lookZ,
//	double upX, double upY, double upZ)
//{
//	renderer->GetActiveCamera()->SetPosition(0, 0, 0);    //相机位置
//	renderer->GetActiveCamera()->SetFocalPoint(lookX, lookY, lookZ);    //焦点位置
//	renderer->GetActiveCamera()->SetViewUp(upX, upY, upZ);    //朝上方向
//	renderer->ResetCamera();
//}
//
//void ViewPositiveX(vtkRenderer* renderer)
//{
//	ViewDirection(renderer, 1, 0, 0, 0, 0, 1);
//}
//
//void ViewNegativeX(vtkRenderer* renderer)
//{
//	ViewDirection(renderer, -1, 0, 0, 0, 0, 1);
//}
//
//void ViewPositiveY(vtkRenderer* renderer)
//{
//	ViewDirection(renderer, 0, 1, 0, 0, 0, 1);
//}
//
//void ViewNegativeY(vtkRenderer* renderer)
//{
//	ViewDirection(renderer, 0, -1, 0, 0, 0, 1);
//}
//
//void ViewPositiveZ(vtkRenderer* renderer)
//{
//	ViewDirection(renderer, 0, 0, 1, 0, 1, 0);
//}
//void ViewNegativeZ(vtkRenderer* renderer)
//{
//	ViewDirection(renderer, 0, 0, -1, 0, 1, 0);
//}
//
//int main(int, char* [])
//{
//	auto namedColors = vtkSmartPointer<vtkNamedColors>::New();
//
//	// X轴上放置立方体
//	auto cubeSource = vtkSmartPointer<vtkCubeSource>::New();
//	cubeSource->SetCenter(1000, 0, 0);
//	cubeSource->SetXLength(1000);
//	cubeSource->SetYLength(1000);
//	cubeSource->SetZLength(1000);
//	cubeSource->Update();
//
//	auto cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
//	auto cubeActor = vtkSmartPointer<vtkActor>::New();
//	cubeActor->SetMapper(cubeMapper);
//	cubeActor->GetProperty()->SetDiffuseColor(
//		namedColors->GetColor3d("Tomato").GetData());
//
//	// Y轴放置圆球
//	auto sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource->SetCenter(0, 1000, 0);
//	sphereSource->SetRadius(500);
//	sphereSource->SetThetaResolution(64);
//	sphereSource->SetPhiResolution(64);
//	sphereSource->Update();
//
//	auto sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
//	auto sphereActor = vtkSmartPointer<vtkActor>::New();
//	sphereActor->SetMapper(sphereMapper);
//	sphereActor->GetProperty()->SetDiffuseColor(
//		namedColors->GetColor3d("Tomato").GetData());
//
//	// Z轴上放置圆锥
//	auto coneSource = vtkSmartPointer<vtkConeSource>::New();
//	coneSource->SetCenter(0, 0, 1000);
//	coneSource->SetRadius(500);
//	coneSource->SetHeight(1000);
//	coneSource->SetDirection(0, 0, 1);
//	coneSource->SetResolution(128);
//	coneSource->Update();
//
//	auto coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	coneMapper->SetInputConnection(coneSource->GetOutputPort());
//	auto coneActor = vtkSmartPointer<vtkActor>::New();
//	coneActor->SetMapper(coneMapper);
//	coneActor->GetProperty()->SetDiffuseColor(
//		namedColors->GetColor3d("Tomato").GetData());
//
//	// Visualize
//	auto renderer = vtkSmartPointer<vtkRenderer>::New();
//	auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->SetSize(500, 500);
//	renderWindow->AddRenderer(renderer);
//
//	auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->SetBackground(namedColors->GetColor3d("SlateGray").GetData());
//	renderer->AddActor(cubeActor);
//	renderer->AddActor(sphereActor);
//	renderer->AddActor(coneActor);
//	renderer->ResetCamera();
//
//	//ViewPositiveX(renderer);
//	//ViewNegativeX(renderer);
//	//ViewPositiveY(renderer);
//	//ViewNegativeY(renderer);
//	//ViewPositiveZ(renderer);
//	ViewNegativeZ(renderer);
//
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//随机平面上色

//#include <vtkVersion.h>
//#include <vtkSmartPointer.h>
//
//#include <vtkActor.h>
//#include <vtkDelaunay2D.h>
//#include <vtkLookupTable.h>
//#include <vtkMath.h>
//#include <vtkPointData.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkVertexGlyphFilter.h>
//#include <vtkXMLPolyDataWriter.h>
//
//// For compatibility with new VTK generic data arrays
//#ifdef vtkGenericDataArray_h
//#define InsertNextTupleValue InsertNextTypedTuple
//#endif
//
//int main(int, char* [])
//{
//	// Create a grid of points (height/terrian map)
//	vtkSmartPointer<vtkPoints> points =
//		vtkSmartPointer<vtkPoints>::New();
//
//	unsigned int GridSize = 20;
//	double xx, yy, zz;
//	//随机生成一个三维坐标，并插入vtkPoints中
//	for (unsigned int x = 0; x < GridSize; x++)
//	{
//		for (unsigned int y = 0; y < GridSize; y++)
//		{
//			xx = x + vtkMath::Random(-.2, .2);
//			yy = y + vtkMath::Random(-.2, .2);
//			zz = vtkMath::Random(-.5, .5);
//			points->InsertNextPoint(xx, yy, zz);
//		}
//	}
//
//	// Add the grid points to a polydata object
//	vtkSmartPointer<vtkPolyData> inputPolyData =
//		vtkSmartPointer<vtkPolyData>::New();
//	inputPolyData->SetPoints(points);
//
//	// Triangulate the grid points
//	//三角化之后，每个面片都是一个三角形，由三个点确定   
//	//一个三角形
//	vtkSmartPointer<vtkDelaunay2D> delaunay =
//		vtkSmartPointer<vtkDelaunay2D>::New();
//
//	//根据vtk的版本来选择数据输入方式
//#if VTK_MAJOR_VERSION <= 5
//	delaunay->SetInput(inputPolyData);
//#else
//	delaunay->SetInputData(inputPolyData);
//#endif
//	delaunay->Update();
//	vtkPolyData* outputPolyData = delaunay->GetOutput();
//
//	//这个数组里面存放了六个坐标范围极值
//	//分别是x轴方向最小最大坐标值
//	//y轴方向最小最大坐标值
//	//Z轴方向最小最大坐标值
//	double bounds[6];
//	outputPolyData->GetBounds(bounds);
//
//	// Find min and max z
//	//数组最后两个值存放的是Z轴方向最小最大坐标值
//	double minz = bounds[4];
//	double maxz = bounds[5];
//
//	std::cout << "minz: " << minz << std::endl;
//	std::cout << "maxz: " << maxz << std::endl;
//
//	// Create the color map
//	vtkSmartPointer<vtkLookupTable> colorLookupTable =
//		vtkSmartPointer<vtkLookupTable>::New();
//	//这里如果不设置大小范围的话，默认的范围是0-1
//	//这里根据自己的代码需要灵活设置颜色映射的范围值
//	colorLookupTable->SetTableRange(minz, maxz);
//	colorLookupTable->Build();
//
//	// Generate the colors for each point based on the color map
//	vtkSmartPointer<vtkUnsignedCharArray> colors =
//		vtkSmartPointer<vtkUnsignedCharArray>::New();
//	colors->SetNumberOfComponents(3);
//	colors->SetName("Colors");
//
//	std::cout << "There are " << outputPolyData->GetNumberOfPoints()
//		<< " points." << std::endl;
//
//	for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
//	{
//		double p[3];
//		outputPolyData->GetPoint(i, p);
//
//		double dcolor[3];
//		//根据z轴的坐标值的数据来获得一个颜色标量值
//		//（使用lookupTable颜色查找表来查找）
//		//查找到的标量值放入dcolor中
//		colorLookupTable->GetColor(p[2], dcolor);
//		std::cout << "dcolor: "
//			<< dcolor[0] << " "
//			<< dcolor[1] << " "
//			<< dcolor[2] << std::endl;
//		unsigned char color[3];
//		//把获取到的颜色标量转换为颜色值，并存入颜色标量colors中
//		for (unsigned int j = 0; j < 3; j++)
//		{
//			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
//		}
//		std::cout << "color: "
//			<< (int)color[0] << " "
//			<< (int)color[1] << " "
//			<< (int)color[2] << std::endl;
//
//		colors->InsertNextTupleValue(color);
//	}
//
//	//给每一个点设置一个颜色，这个颜色是根据z轴的大小来设置的
//	outputPolyData->GetPointData()->SetScalars(colors);
//
//	// Create a mapper and actor
//	vtkSmartPointer<vtkPolyDataMapper> mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//#if VTK_MAJOR_VERSION <= 5
//	mapper->SetInputConnection(outputPolyData->GetProducerPort());
//#else
//	mapper->SetInputData(outputPolyData);
//#endif
//
//	vtkSmartPointer<vtkActor> actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	// Create a renderer, render window, and interactor
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actor to the scene
//	renderer->AddActor(actor);
//	renderer->SetBackground(.1, .2, .3);
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//显示曲率，cxx例子，https://kitware.github.io/vtk-examples/site/Cxx/Visualization/CurvatureBandsWithGlyphs/

//#include <vtkActor.h>
//#include <vtkArrowSource.h>
//#include <vtkBandedPolyDataContourFilter.h>
//#include <vtkCamera.h>
//#include <vtkCleanPolyData.h>
//#include <vtkClipPolyData.h>
//#include <vtkColorSeries.h>
//#include <vtkCurvatures.h>
//#include <vtkElevationFilter.h>
//#include <vtkGlyph3D.h>
//#include <vtkImplicitBoolean.h>
//#include <vtkLookupTable.h>
//#include <vtkMaskPoints.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkParametricFunctionSource.h>
//#include <vtkParametricRandomHills.h>
//#include <vtkParametricTorus.h>
//#include <vtkPlane.h>
//#include <vtkPlaneSource.h>
//#include <vtkPointData.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkReverseSense.h>
//#include <vtkScalarBarActor.h>
//#include <vtkSphereSource.h>
//#include <vtkSuperquadricSource.h>
//#include <vtkTriangleFilter.h>
//#include <vtkVariantArray.h>
//
//#include <algorithm>
//#include <array>
//#include <cctype>
//#include <cmath>
//#include <cstdlib>
//#include <cstring>
//#include <functional>
//#include <iomanip>
//#include <iostream>
//#include <iterator>
//#include <sstream>
//#include <string>
//#include <vector>
//
//enum SURFACE_TYPE
//{
//    TORUS = 0,
//    PARAMETRIC_HILLS,
//    PARAMETRIC_TORUS
//};
//
//namespace {
//    //! Some STL Utilities.
//    class STLHelpers
//    {
//    public:
//        //---------------------------------------------------------------------------
//        STLHelpers()
//        {
//        }
//
//        //---------------------------------------------------------------------------
//        virtual ~STLHelpers()
//        {
//        }
//
//        //-----------------------------------------------------------------------------
//        // Convert a string to lowercase.
//        std::string ToLowercase(const std::string& str)
//        {
//            std::string s;
//            std::transform(str.begin(), str.end(), std::back_inserter(s),
//                (int (*)(int))std::tolower);
//            return s;
//        }
//
//        //-----------------------------------------------------------------------------
//        // Replace all occurrences of old_value in a string with new_value.
//        std::string ReplaceAll(std::string& str, const std::string& old_value,
//            const std::string& new_value)
//        {
//            size_t start_pos = 0;
//            while ((start_pos = str.find(old_value, start_pos)) != std::string::npos)
//            {
//                str.replace(start_pos, old_value.length(), new_value);
//                // It could be that 'new_value' is a substring of 'old_value'.
//                start_pos += new_value.length();
//            }
//            return str;
//        }
//
//        //-----------------------------------------------------------------------------
//        // An implementation of the C++11 next(iter,n) found in the header <iterator>.
//        // ForwardIt must meet the requirements of ForwardIterator.
//        // Return the nth successor of iterator it.
//        template <typename ForwardIt>
//        ForwardIt
//            Next(ForwardIt iter,
//                typename std::iterator_traits<ForwardIt>::difference_type n = 1)
//        {
//            std::advance(iter, n);
//            return iter;
//        }
//
//        //-----------------------------------------------------------------------------
//        // Return true if the iterator points to the last element.
//        template <typename Iter, typename Cont>
//        bool IsLast(Iter iter, const Cont& cont)
//        {
//            return (iter != cont.end()) && (Next(iter) == cont.end());
//        }
//    };
//
//    //-----------------------------------------------------------------------------
//    // Function declarations.
//
//    //! Divide a range into bands
//    /*!
//    @param dR - [min, max] the range that is to be covered by the bands.
//    @param numberOfBands - the number of bands, a positive integer.
//    @param nearestInteger - if True then [floor(min), ceil(max)] is used.
//    @return A List consisting of [min, midpoint, max] for each band.
//    */
//    std::vector<std::vector<double>> MakeBands(double const dR[2],
//        int const& numberOfBands,
//        bool const& nearestInteger);
//
//    //! Divide a range into custom bands
//    /*!
//    You need to specify each band as a list [r1, r2] where r1 < r2 and
//    append these to a list (called x in the implementation).
//    The list should ultimately look
//    like this: x = [[r1, r2], [r2, r3], [r3, r4]...]
//
//    @param dR - [min, max] the range that is to be covered by the bands.
//    @param numberOfBands - the number of bands, a positive integer.
//    @return A List consisting of [min, midpoint, max] for each band.
//    */
//    std::vector<std::vector<double>> MakeCustomBands(double const dR[2],
//        int const& numberOfBands);
//
//    //! Divide a range into integral bands
//    /*!
//    Divide a range into bands
//    @param dR - [min, max] the range that is to be covered by the bands.
//    @return A List consisting of [min, midpoint, max] for each band.
//    */
//    std::vector<std::vector<double>> MakeIntegralBands(double const dR[2]);
//
//    //! Print the bands.
//    /*!
//    @param bands - the bands.
//    */
//    void PrintBands(std::vector<std::vector<double>> const& bands);
//
//    //! Generate elevations over the surface.
//    /*!
//    @param src - the vtkPolyData source.
//    @param elev - the vtkPolyData source with elevations.
//    */
//    void MakeElevations(vtkPolyData* src, vtkPolyData* elev);
//
//    //! Make a torus as the source.
//    /*!
//    @param src - The vtkPolyData source with normal and scalar data.
//    */
//    void MakeTorus(vtkPolyData* src);
//
//    //! Make a parametric torus as the source.
//    /*!
//    @param src - The vtkPolyData source with normal and scalar data.
//    */
//    void MakeParametricTorus(vtkPolyData* src);
//
//    //! Make a parametric hills surface as the source.
//    /*!
//    @param src - The vtkPolyData source with normal and scalar data.
//    */
//    void MakeParametricHills(vtkPolyData* src);
//
//    //! Calculate curvatures.
//    /*!
//    Clip a vtkPolyData source.
//    A cube is made whose size corresponds the the bounds of the source.
//    Then each side is shrunk by the appropriate dx, dy or dz. After
//    this operation the source is clipped by the cube.
//    @param src - the vtkPolyData source
//    @param dx - the amount to clip in the x-direction
//    @param dy - the amount to clip in the y-direction
//    @param dz - the amount to clip in the z-direction
//    @param clipped - clipped vtkPolyData.
//    */
//    void Clipper(vtkPolyData* src, double const& dx, double const& dy,
//        double const& dz, vtkPolyData* clipped);
//
//    //! Calculate curvatures.
//    /*!
//    The source must be triangulated.
//    @param src - the source.
//    @param curv - vtkPolyData with normals and scalar data representing curvatures.
//    */
//    void CalculateCurvatures(vtkPolyData* src, vtkPolyData* curv);
//
//    /*!
//    @param lut - An indexed lookup table.
//    */
//    void MakeLUT(vtkLookupTable* lut);
//
//    //! Create a lookup table with the colors reversed.
//    /*!
//    @param lut - An indexed lookup table.
//    @param lutr - The reversed indexed lookup table.
//    */
//    void ReverseLUT(vtkLookupTable* lut, vtkLookupTable* lutr);
//
//    //! Count the number of scalars in each band.
//    /*!
//    @param bands - the bands.
//    @param src - the vtkPolyData source.
//    @return The frequencies of the scalars in each band.
//    */
//    std::vector<int> Frequencies(std::vector<std::vector<double>> const& bands,
//        vtkPolyData* src);
//
//    //! Print the frequency table.
//    /*!
//    @param freq - the frequencies.
//    */
//    void PrintFrequencies(std::vector<int>& freq);
//
//    //!  Glyph the normals on the surface.
//    /*!
//    @param src - the vtkPolyData source.
//    @param reverseNormals - if True the normals on the surface are reversed.
//    @param glyph - The glyphs.
//    */
//    void MakeGlyphs(vtkPolyData* src, bool const& reverseNormals,
//        vtkGlyph3D* glyph);
//
//    //! Assemble the surface for display.
//    /*!
//    @param st - the surface to display.
//    @param iren - the interactor.
//    */
//    void Display(SURFACE_TYPE st, vtkRenderWindowInteractor* iren);
//    //-----------------------------------------------------------------------------
//
//} // namespace
//
////-----------------------------------------------------------------------------
////! Make and display the surface.
//int main(int, char* [])
//{
//    vtkNew<vtkRenderWindowInteractor> iren;
//    // Select the surface you want displayed.
//    // Display(TORUS, iren);
//    // Display(PARAMETRIC_TORUS, iren);
//    Display(PARAMETRIC_HILLS, iren);
//    iren->Render();
//    iren->Start();
//
//    return EXIT_SUCCESS;
//}
//
//namespace {
//    //-----------------------------------------------------------------------------
//    std::vector<std::vector<double>> MakeBands(double const dR[2],
//        int const& numberOfBands,
//        bool const& nearestInteger)
//    {
//        std::vector<std::vector<double>> bands;
//        if ((dR[1] < dR[0]) || (numberOfBands <= 0))
//        {
//            return bands;
//        }
//        double x[2];
//        for (int i = 0; i < 2; ++i)
//        {
//            x[i] = dR[i];
//        }
//        if (nearestInteger)
//        {
//            x[0] = std::floor(x[0]);
//            x[1] = std::ceil(x[1]);
//        }
//        double dx = (x[1] - x[0]) / static_cast<double>(numberOfBands);
//        std::vector<double> b;
//        b.push_back(x[0]);
//        b.push_back(x[0] + dx / 2.0);
//        b.push_back(x[0] + dx);
//        for (int i = 0; i < numberOfBands; ++i)
//        {
//            bands.push_back(b);
//            for (std::vector<double>::iterator p = b.begin(); p != b.end(); ++p)
//            {
//                *p = *p + dx;
//            }
//        }
//        return bands;
//    }
//
//    //-----------------------------------------------------------------------------
//    std::vector<std::vector<double>> MakeCustomBands(double const dR[2],
//        int const& numberOfBands)
//    {
//        std::vector<std::vector<double>> bands;
//        if ((dR[1] < dR[0]) || (numberOfBands <= 0))
//        {
//            return bands;
//        }
//        // We can do this much better in c++11!
//        double myBands[][2] = { {-0.7, -0.05}, {-0.05, 0},   {0, 0.13},
//                               {0.13, 1.07},  {1.07, 35.4}, {35.4, 37.1} };
//        std::vector<std::vector<double>> x;
//        for (int i = 0; i < 6; ++i)
//        {
//            std::vector<double> tmp(2);
//            tmp[0] = myBands[i][0];
//            tmp[1] = myBands[i][1];
//            x.push_back(tmp);
//        }
//        // Set the minimum to match the range minimum.
//        x[0][0] = dR[0];
//        size_t sz = (static_cast<size_t>(numberOfBands) < x.size())
//            ? static_cast<size_t>(numberOfBands)
//            : x.size();
//        // Adjust the last band.
//        if (x[sz - 1][0] > dR[1])
//        {
//            x[sz - 1][0] = dR[0];
//        }
//        x[sz - 1][1] = dR[1];
//        for (size_t i = 0; i < sz; ++i)
//        {
//            std::vector<double> b(3);
//            b[0] = x[i][0];
//            b[1] = x[i][0] + (x[i][1] - x[i][0]) / 2.0;
//            b[2] = x[i][1];
//            bands.push_back(b);
//        }
//        return bands;
//    }
//
//    //-----------------------------------------------------------------------------
//    std::vector<std::vector<double>> MakeIntegralBands(double const dR[2])
//    {
//        std::vector<std::vector<double>> bands;
//        if (dR[1] < dR[0])
//        {
//            return bands;
//        }
//        double x[2];
//        for (int i = 0; i < 2; ++i)
//        {
//            x[i] = dR[i];
//        }
//        x[0] = std::floor(x[0]);
//        x[1] = std::ceil(x[1]);
//        int numberOfBands = static_cast<int>(std::abs(x[1]) + std::abs(x[0]));
//        return MakeBands(x, numberOfBands, false);
//    }
//
//    //-----------------------------------------------------------------------------
//    void PrintBands(std::vector<std::vector<double>> const& bands)
//    {
//        STLHelpers stlHelpers = STLHelpers();
//        for (std::vector<std::vector<double>>::const_iterator p = bands.begin();
//            p != bands.end(); ++p)
//        {
//            if (p == bands.begin())
//            {
//                std::cout << "[";
//            }
//            for (std::vector<double>::const_iterator q = p->begin(); q != p->end(); ++q)
//            {
//                if (q == p->begin())
//                {
//                    std::cout << "[";
//                }
//                if (!stlHelpers.IsLast(q, *p))
//                {
//                    std::cout << *q << ", ";
//                }
//                else
//                {
//                    std::cout << *q << "]";
//                }
//            }
//            if (!stlHelpers.IsLast(p, bands))
//            {
//                std::cout << ", ";
//            }
//            else
//            {
//                std::cout << "]";
//            }
//        }
//        std::cout << std::endl;
//    }
//
//    //-----------------------------------------------------------------------------
//    void MakeElevations(vtkPolyData* src, vtkPolyData* elev)
//    {
//        double bounds[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//        src->GetBounds(bounds);
//        vtkNew<vtkElevationFilter> elevFilter;
//        elevFilter->SetInputData(src);
//        elevFilter->SetLowPoint(0, bounds[2], 0);
//        elevFilter->SetHighPoint(0, bounds[3], 0);
//        elevFilter->SetScalarRange(bounds[2], bounds[3]);
//        elevFilter->Update();
//        elev->DeepCopy(elevFilter->GetPolyDataOutput());
//    }
//
//    //-----------------------------------------------------------------------------
//    void MakeTorus(vtkPolyData* src)
//    {
//        vtkNew<vtkSuperquadricSource> source;
//        source->SetCenter(0.0, 0.0, 0.0);
//        source->SetCenter(1.0, 1.0, 1.0);
//        source->SetPhiResolution(64);
//        source->SetThetaResolution(64);
//        source->SetThetaRoundness(1);
//        source->SetThickness(0.5);
//        source->SetSize(10);
//        source->SetToroidal(1);
//        source->Update();
//
//        // The quadric is made of strips, so pass it through a triangle filter as
//        // the curvature filter only operates on polys
//        vtkNew<vtkTriangleFilter> tri;
//        tri->SetInputConnection(source->GetOutputPort());
//
//        // The quadric has nasty discontinuities from the way the edges are generated
//        // so let's pass it though a CleanPolyDataFilter and merge any points which
//        // are coincident, or very close
//        vtkNew<vtkCleanPolyData> cleaner;
//        cleaner->SetInputConnection(tri->GetOutputPort());
//        cleaner->SetTolerance(0.005);
//        cleaner->Update();
//
//        vtkNew<vtkPolyData> elev;
//        MakeElevations(cleaner->GetOutput(), elev);
//        CalculateCurvatures(elev, src);
//    }
//
//    //-----------------------------------------------------------------------------
//    void MakeParametricTorus(vtkPolyData* src)
//    {
//        vtkNew<vtkParametricTorus> fn;
//        fn->SetRingRadius(5);
//        fn->SetCrossSectionRadius(2);
//
//        vtkNew<vtkParametricFunctionSource> source;
//        source->SetParametricFunction(fn);
//        source->SetUResolution(50);
//        source->SetVResolution(50);
//        source->SetScalarModeToZ();
//        source->Update();
//        // Name the arrays (not needed in VTK 6.2+ for vtkParametricFunctionSource)
//        source->GetOutput()->GetPointData()->GetNormals()->SetName("Normals");
//        // We have calculated the elevation, just rename the scalars.
//        source->GetOutput()->GetPointData()->GetScalars()->SetName("Elevation");
//        CalculateCurvatures(source->GetOutput(), src);
//    }
//
//    //-----------------------------------------------------------------------------
//    void MakeParametricHills(vtkPolyData* src)
//    {
//        vtkNew<vtkParametricRandomHills> fn;
//        fn->AllowRandomGenerationOn();
//        fn->SetRandomSeed(1);
//        fn->SetNumberOfHills(30);
//        // Make the normals face out of the surface.
//        // Not needed with VTK 8.0 or later.
//        if (strcmp(fn->GetClassName(), "vtkParametricRandomHills") == 0)
//        {
//            fn->ClockwiseOrderingOff();
//        }
//
//        vtkNew<vtkParametricFunctionSource> source;
//        source->SetParametricFunction(fn);
//        source->SetUResolution(50);
//        source->SetVResolution(50);
//        source->SetScalarModeToZ();
//        source->Update();
//        // Name the arrays (not needed in VTK 6.2+ for vtkParametricFunctionSource)
//        source->GetOutput()->GetPointData()->GetNormals()->SetName("Normals");
//        // We have calculated the elevation, just rename the scalars.
//        source->GetOutput()->GetPointData()->GetScalars()->SetName("Elevation");
//        CalculateCurvatures(source->GetOutput(), src);
//    }
//
//    //-----------------------------------------------------------------------------
//    void Clipper(vtkPolyData* src, double const& dx, double const& dy,
//        double const& dz, vtkPolyData* clipped)
//    {
//        double bounds[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//        src->GetBounds(bounds);
//
//        vtkNew<vtkPlane> plane1;
//        plane1->SetOrigin(bounds[0] + dx, 0, 0);
//        plane1->SetNormal(1, 0, 0);
//
//        vtkNew<vtkPlane> plane2;
//        plane2->SetOrigin(bounds[1] - dx, 0, 0);
//        plane2->SetNormal(-1, 0, 0);
//
//        vtkNew<vtkPlane> plane3;
//        plane3->SetOrigin(0, bounds[2] + dy, 0);
//        plane3->SetNormal(0, 1, 0);
//
//        vtkNew<vtkPlane> plane4;
//        plane4->SetOrigin(0, bounds[3] - dy, 0);
//        plane4->SetNormal(0, -1, 0);
//
//        vtkNew<vtkPlane> plane5;
//        plane5->SetOrigin(0, 0, bounds[4] + dz);
//        plane5->SetNormal(0, 0, 1);
//
//        vtkNew<vtkPlane> plane6;
//        plane6->SetOrigin(0, 0, bounds[5] - dz);
//        plane6->SetNormal(0, 0, -1);
//
//        vtkNew<vtkImplicitBoolean> clipFunction;
//        clipFunction->SetOperationTypeToUnion();
//        clipFunction->AddFunction(plane1);
//        clipFunction->AddFunction(plane2);
//        clipFunction->AddFunction(plane3);
//        clipFunction->AddFunction(plane4);
//        clipFunction->AddFunction(plane5);
//        clipFunction->AddFunction(plane6);
//
//        // Create clipper for the random hills
//        vtkNew<vtkClipPolyData> clipper;
//        clipper->SetClipFunction(clipFunction);
//        clipper->SetInputData(src);
//        clipper->GenerateClipScalarsOff();
//        clipper->GenerateClippedOutputOff();
//        // clipper->GenerateClippedOutputOn();
//        clipper->Update();
//        clipped->DeepCopy(clipper->GetOutput());
//    }
//
//    //-----------------------------------------------------------------------------
//    void CalculateCurvatures(vtkPolyData* src, vtkPolyData* curv)
//    {
//        // Calculate the curvature.
//        vtkNew<vtkCurvatures> curvature;
//        curvature->SetCurvatureTypeToGaussian();
//        curvature->SetInputData(src);
//        curvature->Update();
//        curv->DeepCopy(curvature->GetOutput());
//    }
//
//    //-----------------------------------------------------------------------------
//    void MakeLUT(vtkLookupTable* lut)
//    {
//        // Make the lookup table.
//        vtkNew<vtkColorSeries> colorSeries;
//        // Select a color scheme.
//        int colorSeriesEnum;
//        // colorSeriesEnum = colorSeries->BREWER_DIVERGING_BROWN_BLUE_GREEN_9;
//        // colorSeriesEnum = colorSeries->BREWER_DIVERGING_SPECTRAL_10;
//        // colorSeriesEnum = colorSeries->BREWER_DIVERGING_SPECTRAL_3;
//        // colorSeriesEnum = colorSeries->BREWER_DIVERGING_PURPLE_ORANGE_9;
//        // colorSeriesEnum = colorSeries->BREWER_SEQUENTIAL_BLUE_PURPLE_9;
//        // colorSeriesEnum = colorSeries->BREWER_SEQUENTIAL_BLUE_GREEN_9;
//        colorSeriesEnum = colorSeries->BREWER_QUALITATIVE_SET3;
//        // colorSeriesEnum = colorSeries->CITRUS;
//        colorSeries->SetColorScheme(colorSeriesEnum);
//        colorSeries->BuildLookupTable(lut);
//        lut->SetNanColor(0, 0, 0, 1);
//    }
//
//    //! Create a lookup table with the colors reversed.
//    /*!
//    @param lut - An indexed lookup table.
//    @param lutr - The reversed indexed lookup table.
//    */
//    //-----------------------------------------------------------------------------
//    void ReverseLUT(vtkLookupTable* lut, vtkLookupTable* lutr)
//    {
//        // First do a deep copy just to get the whole structure
//        // and then reverse the colors and annotations.
//        lutr->DeepCopy(lut);
//        vtkIdType t = lut->GetNumberOfTableValues() - 1;
//        for (vtkIdType i = t; i >= 0; --i)
//        {
//            double rgba[4] = { 0.0, 0.0, 0.0, 0.0 };
//            lut->GetColor(i, rgba);
//            rgba[3] = lut->GetOpacity(i);
//            lutr->SetTableValue(t - i, rgba);
//        }
//        t = lut->GetNumberOfAnnotatedValues() - 1;
//        for (vtkIdType i = t; i >= 0; --i)
//        {
//            lutr->SetAnnotation(t - i, lut->GetAnnotation(i));
//        }
//    }
//
//    //-----------------------------------------------------------------------------
//    std::vector<int> Frequencies(std::vector<std::vector<double>> const& bands,
//        vtkPolyData* src)
//    {
//        std::vector<int> freq(bands.size(), 0);
//        vtkIdType tuples = src->GetPointData()->GetScalars()->GetNumberOfTuples();
//        for (int i = 0; i < tuples; ++i)
//        {
//            double* x = src->GetPointData()->GetScalars()->GetTuple(i);
//            for (size_t j = 0; j < bands.size(); ++j)
//            {
//                if (*x <= bands[j][2])
//                {
//                    freq[j] = freq[j] + 1;
//                    break;
//                }
//            }
//        }
//        return freq;
//    }
//
//    //-----------------------------------------------------------------------------
//    void PrintFrequencies(std::vector<int>& freq)
//    {
//        STLHelpers stlHelpers = STLHelpers();
//        int i = 0;
//        for (std::vector<int>::const_iterator p = freq.begin(); p != freq.end(); ++p)
//        {
//            if (p == freq.begin())
//            {
//                std::cout << "[";
//            }
//            if (stlHelpers.IsLast(p, freq))
//            {
//                std::cout << i << ": " << *p << "]";
//            }
//            else
//            {
//                std::cout << i << ": " << *p << ", ";
//            }
//            ++i;
//        }
//        std::cout << endl;
//    }
//
//    //-----------------------------------------------------------------------------
//    void MakeGlyphs(vtkPolyData* src, bool const& reverseNormals, vtkGlyph3D* glyph)
//    {
//        // Sometimes the contouring algorithm can create a volume whose gradient
//        // vector and ordering of polygon(using the right hand rule) are
//        // inconsistent. vtkReverseSense cures this problem.
//        vtkNew<vtkReverseSense> reverse;
//        vtkNew<vtkMaskPoints> maskPts;
//        maskPts->SetOnRatio(5);
//        maskPts->RandomModeOn();
//        if (reverseNormals)
//        {
//            reverse->SetInputData(src);
//            reverse->ReverseCellsOn();
//            reverse->ReverseNormalsOn();
//            maskPts->SetInputConnection(reverse->GetOutputPort());
//        }
//        else
//        {
//            maskPts->SetInputData(src);
//        }
//
//        // Source for the glyph filter
//        vtkNew<vtkArrowSource> arrow;
//        arrow->SetTipResolution(16);
//        arrow->SetTipLength(0.3);
//        arrow->SetTipRadius(0.1);
//
//        glyph->SetSourceConnection(arrow->GetOutputPort());
//        glyph->SetInputConnection(maskPts->GetOutputPort());
//        glyph->SetVectorModeToUseNormal();
//        glyph->SetScaleFactor(1);
//        glyph->SetColorModeToColorByVector();
//        glyph->SetScaleModeToScaleByVector();
//        glyph->OrientOn();
//        glyph->Update();
//    }
//
//    //-----------------------------------------------------------------------------
//    void Display(SURFACE_TYPE st, vtkRenderWindowInteractor* iren)
//    {
//
//        vtkNew<vtkNamedColors> colors;
//
//        // Set the background color.
//        std::array<unsigned char, 4> bkg{ {179, 204, 255, 255} };
//        colors->SetColor("BkgColor", bkg.data());
//
//        // ------------------------------------------------------------
//        // Create the surface, lookup tables, contour filter etc.
//        // ------------------------------------------------------------
//        vtkNew<vtkPolyData> src;
//        switch (st)
//        {
//        case TORUS: {
//            MakeTorus(src);
//            break;
//        }
//        case PARAMETRIC_TORUS: {
//            MakeParametricTorus(src);
//            break;
//        }
//        case PARAMETRIC_HILLS: {
//            vtkNew<vtkPolyData> hills;
//            MakeParametricHills(hills);
//            Clipper(hills, 0.5, 0.5, 0.0, src);
//            break;
//        }
//        default: {
//            std::cout << "No surface specified." << std::endl;
//            return;
//        }
//        }
//        //  Here we are assuming that the active scalars are the curvatures.
//        // in the parametric surfaces, so change the name.
//        char* curvatureName = src->GetPointData()->GetScalars()->GetName();
//        //  Use this range to color the glyphs for the normals by elevation.
//        src->GetPointData()->SetActiveScalars("Elevation");
//        double scalarRangeElevation[2];
//        src->GetScalarRange(scalarRangeElevation);
//        src->GetPointData()->SetActiveScalars(curvatureName);
//        double scalarRangeCurvatures[2];
//        src->GetScalarRange(scalarRangeCurvatures);
//        double scalarRange[2];
//        scalarRange[0] = scalarRangeCurvatures[0];
//        scalarRange[1] = scalarRangeCurvatures[1];
//
//        vtkNew<vtkLookupTable> lut;
//        MakeLUT(lut);
//        vtkIdType numberOfBands = lut->GetNumberOfTableValues();
//        std::vector<std::vector<double>> bands;
//        if (st == PARAMETRIC_HILLS)
//        {
//            // Comment this out if you want to see how allocating
//            // equally spaced bands works.
//            bands = MakeCustomBands(scalarRange, numberOfBands);
//            // Adjust the number of table values
//            numberOfBands = static_cast<vtkIdType>(bands.size());
//            lut->SetNumberOfTableValues(numberOfBands);
//        }
//        else
//        {
//            bands = MakeBands(scalarRange, numberOfBands, false);
//        }
//        lut->SetTableRange(scalarRange);
//
//        // PrintBands(bands);
//
//        // Let's do a frequency table.
//        // The number of scalars in each band.
//        // std::vector<int> freq = Frequencies(bands, src);
//        // PrintFrequencies(freq);
//
//        // We will use the midpoint of the band as the label.
//        std::vector<std::string> labels;
//        for (std::vector<std::vector<double>>::const_iterator p = bands.begin();
//            p != bands.end(); ++p)
//        {
//            std::ostringstream os;
//            os << std::fixed << std::setw(6) << std::setprecision(2) << (*p)[1];
//            labels.push_back(os.str());
//        }
//
//        // Annotate
//        vtkNew<vtkVariantArray> values;
//        for (size_t i = 0; i < labels.size(); ++i)
//        {
//            values->InsertNextValue(vtkVariant(labels[i]));
//        }
//        for (vtkIdType i = 0; i < values->GetNumberOfTuples(); ++i)
//        {
//            lut->SetAnnotation(i, values->GetValue(i).ToString());
//        }
//
//        // Create a lookup table with the colors reversed.
//        vtkNew<vtkLookupTable> lutr;
//        ReverseLUT(lut, lutr);
//
//        // Create the contour bands.
//        vtkNew<vtkBandedPolyDataContourFilter> bcf;
//        bcf->SetInputData(src);
//        // Use either the minimum or maximum value for each band.
//        int i = 0;
//        for (std::vector<std::vector<double>>::const_iterator p = bands.begin();
//            p != bands.end(); ++p)
//        {
//            bcf->SetValue(i, (*p)[2]);
//            ++i;
//        }
//        // We will use an indexed lookup table.
//        bcf->SetScalarModeToIndex();
//        bcf->GenerateContourEdgesOn();
//
//        // Generate the glyphs on the original surface.
//        vtkNew<vtkGlyph3D> glyph;
//        MakeGlyphs(src, false, glyph);
//
//        // ------------------------------------------------------------
//        // Create the mappers and actors
//        // ------------------------------------------------------------
//
//        vtkNew<vtkPolyDataMapper> srcMapper;
//        srcMapper->SetInputConnection(bcf->GetOutputPort());
//        srcMapper->SetScalarRange(scalarRange);
//        srcMapper->SetLookupTable(lut);
//        srcMapper->SetScalarModeToUseCellData();
//
//        vtkNew<vtkActor> srcActor;
//        srcActor->SetMapper(srcMapper);
//        srcActor->RotateX(-45);
//        srcActor->RotateZ(45);
//
//        // Create contour edges
//        vtkNew<vtkPolyDataMapper> edgeMapper;
//        edgeMapper->SetInputData(bcf->GetContourEdgesOutput());
//        edgeMapper->SetResolveCoincidentTopologyToPolygonOffset();
//
//        vtkNew<vtkActor> edgeActor;
//        edgeActor->SetMapper(edgeMapper);
//        edgeActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
//        edgeActor->RotateX(-45);
//        edgeActor->RotateZ(45);
//
//        vtkNew<vtkPolyDataMapper> glyphMapper;
//        glyphMapper->SetInputConnection(glyph->GetOutputPort());
//        glyphMapper->SetScalarModeToUsePointFieldData();
//        glyphMapper->SetColorModeToMapScalars();
//        glyphMapper->ScalarVisibilityOn();
//        glyphMapper->SelectColorArray("Elevation");
//        // Colour by scalars.
//        // The default lookup table is used but you can
//        // use whatever lookup table you like.
//        glyphMapper->SetScalarRange(scalarRangeElevation);
//
//        vtkNew<vtkActor> glyphActor;
//        glyphActor->SetMapper(glyphMapper);
//        glyphActor->RotateX(-45);
//        glyphActor->RotateZ(45);
//
//        // Add a scalar bar->
//        vtkNew<vtkScalarBarActor> scalarBar;
//        // This LUT puts the lowest value at the top of the scalar bar.
//        // scalarBar->SetLookupTable(lut);
//        // Use this LUT if you want the highest value at the top.
//        scalarBar->SetLookupTable(lutr);
//        scalarBar->SetTitle("Gaussian\nCurvature");
//
//        // ------------------------------------------------------------
//        // Create the RenderWindow, Renderer and Interactor
//        // ------------------------------------------------------------
//        vtkNew<vtkRenderer> ren;
//        vtkNew<vtkRenderWindow> renWin;
//
//        renWin->AddRenderer(ren);
//        iren->SetRenderWindow(renWin);
//
//        // add actors
//        ren->AddViewProp(srcActor);
//        ren->AddViewProp(edgeActor);
//        ren->AddViewProp(glyphActor);
//        ren->AddActor2D(scalarBar);
//
//        ren->SetBackground(colors->GetColor3d("BkgColor").GetData());
//        renWin->SetSize(800, 800);
//        renWin->SetWindowName("CurvatureBandsWithGlyphs");
//        renWin->Render();
//
//        ren->GetActiveCamera()->Zoom(1.5);
//    }
//
//} // namespace



//显示曲率和坐标色卡，装配.stl
/**********************************************************************

Copyright (c) Mr.Bin. All rights reserved.
For more information visit: http://blog.csdn.net/webzhuce

**********************************************************************/
/*
#include <vtkSmartPointer.h>
#include <vtkCurvatures.h>
#include <vtkPolyDataReader.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkColorSeries.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include<vtkUnstructuredGridReader.h>
#include<vtkSTLReader.h>

//test data : ../data/fran_cut.vtk
int main(int argc, char* argv[])
{
	//vtkNew<vtkPolyDataReader> reader;
	//reader->SetFileName("zhuangpei.stl");
	//reader->Update();
	vtkSmartPointer<vtkPolyData> input1 = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName("luntai.stl");
	reader->Update();
	input1->DeepCopy(reader->GetOutput());

	vtkNew<vtkCurvatures> curvaturesfilter;
	curvaturesfilter->SetInputConnection(reader->GetOutputPort());
	curvaturesfilter->SetCurvatureTypeToMinimum();
	//curvaturesfilter->SetCurvatureTypeToMaximum();
	//curvaturesfilter->SetCurvatureTypeToGaussian();
	//curvaturesfilter->SetCurvatureTypeToMean();
	curvaturesfilter->Update();

	double scalarrange[2];
	curvaturesfilter->GetOutput()->GetScalarRange(scalarrange);
	vtkNew<vtkLookupTable> lut;
	lut->SetHueRange(0.0, 0.6);
	lut->SetAlphaRange(1.0, 1.0);
	lut->SetValueRange(1.0, 1.0);
	lut->SetSaturationRange(1.0, 1.0);
	lut->SetNumberOfTableValues(256);
	lut->SetRange(scalarrange);
	lut->Build();

	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(curvaturesfilter->GetOutputPort());
	mapper->SetLookupTable(lut);
	mapper->SetScalarRange(scalarrange);
	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	vtkNew<vtkScalarBarActor> scalarbar;
	scalarbar->SetLookupTable(mapper->GetLookupTable());
	scalarbar->SetTitle(curvaturesfilter->GetOutput()->GetPointData()->GetScalars()->GetName());
	scalarbar->SetNumberOfLabels(5);
	vtkNew<vtkRenderer> renderer;
	renderer->AddActor(actor);
	renderer->AddActor2D(scalarbar);
	renderer->SetBackground(1.0, 1.0, 1.0);
	vtkNew<vtkRenderWindow> renderwindow;
	renderwindow->AddRenderer(renderer);
	renderwindow->SetSize(640, 480);
	renderwindow->Render();
	renderwindow->SetWindowName("PolyDataCurvature");
	vtkNew<vtkRenderWindowInteractor> renderwindowinteractor;
	renderwindowinteractor->SetRenderWindow(renderwindow);
	renderwindow->Render();
	renderwindowinteractor->Start();

	return EXIT_SUCCESS;
}*/


//创建趋势箭头
/*
#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkAssignAttribute.h>
#include <vtkCamera.h>
#include <vtkExtractEdges.h>
#include <vtkGlyph3D.h>
#include <vtkGradientFilter.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTubeFilter.h>
#include <vtkUnstructuredGridReader.h>

int main(int argc, char* argv[])
{
	// Create the reader for the data.
	// This is the data that will be rendered.
	std::string filename = "Test.vtk";
	std::cout << "Loading " << filename.c_str() << std::endl;
	vtkNew<vtkUnstructuredGridReader> reader;
	reader->SetFileName(filename.c_str());

	vtkNew<vtkExtractEdges> edges;
	edges->SetInputConnection(reader->GetOutputPort());

	vtkNew<vtkTubeFilter> tubes;//管道
	tubes->SetInputConnection(edges->GetOutputPort());
	tubes->SetRadius(0.0625);
	tubes->SetVaryRadiusToVaryRadiusOff();
	tubes->SetNumberOfSides(32);

	vtkNew<vtkPolyDataMapper> tubesMapper;
	tubesMapper->SetInputConnection(tubes->GetOutputPort());
	tubesMapper->SetScalarRange(0.0, 26.0);

	vtkNew<vtkActor> tubesActor;
	tubesActor->SetMapper(tubesMapper);

	vtkNew<vtkGradientFilter> gradients;//梯度
	gradients->SetInputConnection(reader->GetOutputPort());

	vtkNew<vtkAssignAttribute> vectors;
	vectors->SetInputConnection(gradients->GetOutputPort());
	vectors->Assign("Gradients", vtkDataSetAttributes::VECTORS,
		vtkAssignAttribute::POINT_DATA);

	vtkNew<vtkArrowSource> arrow;

	vtkNew<vtkGlyph3D> glyphs;//箭头
	glyphs->SetInputConnection(0, vectors->GetOutputPort());
	glyphs->SetInputConnection(1, arrow->GetOutputPort());
	glyphs->ScalingOn();
	glyphs->SetScaleModeToScaleByVector();
	glyphs->SetScaleFactor(0.25);
	glyphs->OrientOn();
	glyphs->ClampingOff();
	glyphs->SetVectorModeToUseVector();
	glyphs->SetIndexModeToOff();

	vtkNew<vtkPolyDataMapper> glyphMapper;
	glyphMapper->SetInputConnection(glyphs->GetOutputPort());
	glyphMapper->ScalarVisibilityOff();

	vtkNew<vtkActor> glyphActor;
	glyphActor->SetMapper(glyphMapper);

	vtkNew<vtkRenderer> renderer;
	renderer->AddActor(tubesActor);
	renderer->AddActor(glyphActor);
	renderer->SetBackground(0.328125, 0.347656, 0.425781);

	vtkNew<vtkRenderWindow> renwin;
	renwin->AddRenderer(renderer);
	renwin->SetSize(350, 500);
	renwin->SetWindowName("GradientFilter");

	renderer->ResetCamera();
	vtkCamera* camera = renderer->GetActiveCamera();
	camera->Elevation(-80.0);
	camera->OrthogonalizeViewUp();
	camera->Azimuth(135.0);

	vtkNew<vtkRenderWindowInteractor> iren;
	iren->SetRenderWindow(renwin);
	renwin->Render();
	iren->Initialize();
	iren->Start();

	return EXIT_SUCCESS;
}
*/



//读.vtk文件，染色，显示色卡，选择显示标量

//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkLookupTable.h>
//#include <algorithm>
//#include <array>
//#include <string>
//#include <vtkScalarBarActor.h>
//#include <vtkAxesActor.h>
//#include <vtkOrientationMarkerWidget.h>
//#include <vtkDataSetAttributes.h>
//#include <vtkPolyData.h>
//#include <vtkPointData.h>
//#include <vtkCellData.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	vtkNew<vtkInteractorStyleTrackballCamera>style;
//	interactor->SetInteractorStyle(style);
//	interactor->SetRenderWindow(renderWindow);
//
//	renderer->SetBackground(colors->GetColor3d("White").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	reader->ReadAllScalarsOn();//获取所有的标量数据
//	reader->ReadAllVectorsOn();
//	reader->ReadAllNormalsOn();
//	reader->ReadAllTensorsOn();
//	reader->ReadAllColorScalarsOn();
//	reader->ReadAllTCoordsOn();
//	reader->ReadAllFieldsOn();
//
//	reader->SetFileName("Test04.vtk");
//	reader->GetOutput()->Register(reader);
//	reader->Update();
//
//
//	int nNumScalar = reader->GetNumberOfScalarsInFile();//获取标量类型数
//	cout << nNumScalar << endl;
//
//
//	//std::cout << "Loading: " << argv[1] << std::endl;
//	//auto unstructuredGrid = ReadUnstructuredGrid(std::string(argv[1]));
//	vtkSmartPointer<vtkPolyData> UnstructuredGrid =
//		vtkSmartPointer<vtkPolyData>::New();
//	auto unstructuredGrid = reader->GetOutput();
//	for (int i = 0; i < nNumScalar; i++)
//	{
//		cout << reader->GetScalarsNameInFile(i) << endl;
//		//cout << reader->GetScalarsNameInFile(1) << endl;
//		//cout << reader->GetScalarsNameInFile(2) << endl;
//		//cout << reader->GetScalarsNameInFile(3) << endl;
//		//cout << reader->GetScalarsNameInFile(4) << endl;
//	}
//	int index = 0;
//
//	//reader->GetOutput()->GetPointData()->SetActiveScalars(reader->GetScalarsNameInFile(index));//设置标量名称，即渲染哪个标量
//	reader->GetOutput()->GetCellData()->SetActiveScalars(reader->GetScalarsNameInFile(index));
//
//
//	vtkNew<vtkLookupTable> lut1;
//	lut1->SetHueRange(.667, 0);// 设定HSV颜色范围，色调H取值范围为0°～360°，从红色开始按逆时针方向计算，红色为0°/0.0，绿色为120°/0.34,蓝色为240°/0.67
//	// Visualize
//	vtkNew<vtkDataSetMapper> mapper;
//	mapper->SetInputData(unstructuredGrid);
//	//mapper->ScalarVisibilityOff();
//	mapper->SetScalarRange(unstructuredGrid->GetScalarRange());
//	mapper->SetLookupTable(lut1);
//	mapper->SetColorModeToMapScalars();
//	//mapper->SetScalarModeToUseCellData();
//	mapper->SetScalarModeToDefault();
//
//	cout << unstructuredGrid->GetScalarRange()[0] << endl;
//	cout << unstructuredGrid->GetScalarRange()[1] << endl;
//
//	vtkNew<vtkScalarBarActor> scalarbar;
//	scalarbar->SetLookupTable(mapper->GetLookupTable());
//	//scalarbar->SetTitle(curvaturesfilter->GetOutput()->GetPointData()->GetScalars()->GetName());
//	scalarbar->SetNumberOfLabels(5);
//	scalarbar->SetTitle(reader->GetScalarsNameInFile(index));
//	renderer->AddActor2D(scalarbar);
//
//
//	vtkNew<vtkAxesActor> axes;
//
//	vtkNew<vtkOrientationMarkerWidget> widget;
//	double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
//	colors->GetColor("Carrot", rgba);
//	widget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
//	widget->SetOrientationMarker(axes);
//	widget->SetInteractor(interactor);
//	widget->SetViewport(0.0, 0.0, 0.4, 0.4);
//	widget->SetEnabled(1);
//	widget->InteractiveOn();
//
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->GetProperty()->EdgeVisibilityOn();//显示网格
//	//actor->GetProperty()->SetSpecular(.3);
//	//actor->GetProperty()->SetSpecularPower(30);
//	//actor->GetProperty()->EdgeVisibilityOn();
//	actor->GetProperty()->SetOpacity(.5);
//
//	renderer->AddActor(actor);
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}


//2021/4/16 忘了是什么了
//#include <vtkVersion.h>
//#include <vtkSmartPointer.h>
//#include <vtkPointPicker.h>
//#include <vtkSphereSource.h>
//#include <vtkGlyph3D.h>
//#include <vtkPointData.h>
//#include <vtkIdTypeArray.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkRendererCollection.h>
//#include <vtkProperty.h>
//#include <vtkPlanes.h>
//#include <vtkObjectFactory.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkPolyData.h>
//#include <vtkPointSource.h>
//#include <vtkInteractorStyleTrackballActor.h>
//#include <vtkAreaPicker.h>
//#include <vtkExtractGeometry.h>
//#include <vtkDataSetMapper.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkVertexGlyphFilter.h>
//#include <vtkIdFilter.h>
//
//// Define interaction style
//class InteractorStyle2 : public vtkInteractorStyleTrackballActor
//{
//public:
//	static InteractorStyle2* New();
//	vtkTypeMacro(InteractorStyle2, vtkInteractorStyleTrackballActor);
//
//	InteractorStyle2()
//	{
//		this->Move = false;
//		this->PointPicker = vtkSmartPointer<vtkPointPicker>::New();
//
//		// Setup ghost glyph
//		vtkSmartPointer<vtkPoints> points =
//			vtkSmartPointer<vtkPoints>::New();
//		points->InsertNextPoint(0, 0, 0);
//		this->MovePolyData = vtkSmartPointer<vtkPolyData>::New();
//		this->MovePolyData->SetPoints(points);
//		this->MoveGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
//#if VTK_MAJOR_VERSION <= 5
//		this->MoveGlyphFilter->SetInputConnection(
//			this->MovePolyData->GetProducerPort());
//#else
//		this->MoveGlyphFilter->SetInputData(this->MovePolyData);
//#endif
//		this->MoveGlyphFilter->Update();
//
//		this->MoveMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//		this->MoveMapper->SetInputConnection(this->MoveGlyphFilter->GetOutputPort());
//
//		this->MoveActor = vtkSmartPointer<vtkActor>::New();
//		this->MoveActor->SetMapper(this->MoveMapper);
//		this->MoveActor->VisibilityOff();
//		this->MoveActor->GetProperty()->SetPointSize(10);
//		this->MoveActor->GetProperty()->SetColor(1, 0, 0);
//	}
//
//	void OnMouseMove()
//	{
//		if (!this->Move)
//		{
//			return;
//		}
//
//		vtkInteractorStyleTrackballActor::OnMouseMove();
//
//	}
//
//	void OnMiddleButtonUp()
//	{
//		this->EndPan();
//
//		this->Move = false;
//		this->MoveActor->VisibilityOff();
//
//		this->Data->GetPoints()->SetPoint(this->SelectedPoint, this->MoveActor->GetPosition());
//		this->Data->Modified();
//		this->GetCurrentRenderer()->Render();
//		this->GetCurrentRenderer()->GetRenderWindow()->Render();
//
//	}
//	void OnMiddleButtonDown()
//	{
//		// Get the selected point
//		int x = this->Interactor->GetEventPosition()[0];
//		int y = this->Interactor->GetEventPosition()[1];
//		this->FindPokedRenderer(x, y);
//
//		this->PointPicker->Pick(this->Interactor->GetEventPosition()[0],
//			this->Interactor->GetEventPosition()[1],
//			0,  // always zero.
//			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
//
//		if (this->PointPicker->GetPointId() >= 0)
//		{
//			this->StartPan();
//			this->MoveActor->VisibilityOn();
//			this->Move = true;
//			this->SelectedPoint = this->PointPicker->GetPointId();
//
//			std::cout << "Dragging point " << this->SelectedPoint << std::endl;
//
//			double p[3];
//			this->Data->GetPoint(this->SelectedPoint, p);
//			std::cout << "p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//			this->MoveActor->SetPosition(p);
//
//			this->GetCurrentRenderer()->AddActor(this->MoveActor);
//			this->InteractionProp = this->MoveActor;
//		}
//	}
//
//	vtkPolyData* Data;
//	vtkPolyData* GlyphData;
//
//	vtkSmartPointer<vtkPolyDataMapper> MoveMapper;
//	vtkSmartPointer<vtkActor> MoveActor;
//	vtkSmartPointer<vtkPolyData> MovePolyData;
//	vtkSmartPointer<vtkVertexGlyphFilter> MoveGlyphFilter;
//
//	vtkSmartPointer<vtkPointPicker> PointPicker;
//
//	bool Move;
//	vtkIdType SelectedPoint;
//};
//vtkStandardNewMacro(InteractorStyle2);
//
//int main(int, char* [])
//{
//	vtkSmartPointer<vtkPoints> points =
//		vtkSmartPointer<vtkPoints>::New();
//	points->InsertNextPoint(0, 0, 0);
//	points->InsertNextPoint(1, 0, 0);
//	points->InsertNextPoint(2, 0, 0);
//
//	vtkSmartPointer<vtkPolyData> input =
//		vtkSmartPointer<vtkPolyData>::New();
//	input->SetPoints(points);
//
//	vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//#if VTK_MAJOR_VERSION <= 5
//	glyphFilter->SetInputConnection(input->GetProducerPort());
//#else
//	glyphFilter->SetInputData(input);
//#endif
//	glyphFilter->Update();
//
//	// Create a mapper and actor
//	vtkSmartPointer<vtkPolyDataMapper> mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(glyphFilter->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetPointSize(10);
//
//	// Visualize
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(actor);
//	//renderer->SetBackground(1,1,1); // Background color white
//
//	renderWindow->Render();
//
//	vtkSmartPointer<InteractorStyle2> style =
//		vtkSmartPointer<InteractorStyle2>::New();
//	renderWindowInteractor->SetInteractorStyle(style);
//	style->Data = input;
//
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}



//显示坐标轴，但是会随着视角旋转改变坐标轴位置
/*
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>

int main(int, char* [])
{
	vtkNew<vtkNamedColors> colors;

	//vtkNew<vtkSphereSource> sphereSource;
	//sphereSource->SetCenter(0.0, 0.0, 0.0);
	//sphereSource->SetRadius(0.5);

	// create a mapper
	//vtkNew<vtkPolyDataMapper> sphereMapper;
	//sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

	// create an actor
	//vtkNew<vtkActor> sphereActor;
	//sphereActor->SetMapper(sphereMapper);

	// a renderer and render window
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetWindowName("Axes");
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(300, 300);

	// an interactor
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// add the actors to the scene
	//renderer->AddActor(sphereActor);
	renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

	vtkNew<vtkTransform> transform;
	transform->Translate(0.0, 0.0, 0.0);

	vtkNew<vtkAxesActor> axes;

	// The axes are positioned with a user transform
	//axes->SetUserTransform(transform);

	// properties of the axes labels can be set as follows
	// this sets the x axis label to red
	// axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(
	//   colors->GetColor3d("Red").GetData());

	// the actual text of the axis label can be changed:
	// axes->SetXAxisLabelText("test");

	renderer->AddActor(axes);

	renderer->GetActiveCamera()->Azimuth(50);
	renderer->GetActiveCamera()->Elevation(-30);

	renderer->ResetCamera();
	renderWindow->SetWindowName("Axes");
	renderWindow->Render();

	// begin mouse interaction
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}*/



//单独坐标轴
/*
#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPropAssembly.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>

int main(int, char* [])
{
	vtkNew<vtkNamedColors> colors;

	//vtkNew<vtkSphereSource> sphereSource;
	//sphereSource->SetCenter(0.0, 0.0, 0.0);
	//sphereSource->SetRadius(1.0);
	//sphereSource->Update();

	//vtkPolyData* polydata = sphereSource->GetOutput();

	// Create a mapper
	//vtkNew<vtkPolyDataMapper> mapper;
	//mapper->SetInputData(polydata);

	// Create an actor
	//vtkNew<vtkActor> actor;
	//actor->SetMapper(mapper);
	//actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

	// A renderer and render window
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetWindowName("DisplayCoordinateAxes");
	renderWindow->AddRenderer(renderer);

	// An interactor
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	//renderer->AddActor(actor);
	//renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

	vtkNew<vtkAxesActor> axes;

	vtkNew<vtkOrientationMarkerWidget> widget;
	double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
	colors->GetColor("Carrot", rgba);
	widget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
	widget->SetOrientationMarker(axes);
	widget->SetInteractor(renderWindowInteractor);
	widget->SetViewport(0.0, 0.0, 0.4, 0.4);
	widget->SetEnabled(1);
	widget->InteractiveOn();

	renderer->GetActiveCamera()->Azimuth(50);
	renderer->GetActiveCamera()->Elevation(-30);

	renderer->ResetCamera();
	renderWindow->Render();

	// Begin mouse interaction
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}*/


//画直线，正常的鼠标缩放、移动等视图控制
/*
#include <vtkSmartPointer.h>
#include <vtkAnimationCue.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkCommand.h>
#include <vtkAnimationScene.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLineSource.h>
#include "vtkInteractorStyleTrackballCamera.h"
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>

class vtkCustomAnimationCue : public vtkAnimationCue
{
public:
	static vtkCustomAnimationCue* New()
	{
		vtkCustomAnimationCue* p = new vtkCustomAnimationCue();
		return p;
	}
	vtkRenderWindow* RenWin;  //活跃的渲染器
	vtkLineSource* line;
	vtkPoints* points;
protected:
	double ptemp[3]; //临时点的储存
	vtkCustomAnimationCue()
	{
		this->RenWin = 0;
		this->line = 0;
		this->points = 0;
		//对这个临时点的初始位置进行初始化
		ptemp[0] = 0;//X坐标
		ptemp[1] = 1;//Y坐标
		ptemp[2] = 10; //Z坐标
	}
	virtual void TickInternal(double currenttime, double deltatime, double clocktime)  //这个函数即在每一帧调用的重新绘制函数
	{
		ptemp[2] += 0.1;//对于每一帧，把初始值的y坐标加0.1
		this->points->InsertNextPoint(ptemp);
		this->line->Modified();//放入点集中并对直线的数据进行修改
		this->RenWin->Render();//重新渲染
	}
};
int main()
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
	//生成四个点，得到三段连起来的直线
	double p0[3] = { 0, 0.0, 7.0 };
	double p1[3] = { 3.0, 10.0, 7.0 };
	double p2[3] = { 63.0,15.0,7.0 };
	double p3[3] = { 0,18,10 };
	points->InsertNextPoint(p0);
	points->InsertNextPoint(p1);
	points->InsertNextPoint(p2);
	points->InsertNextPoint(p3);
	line->SetPoints(points);
	line->Update();


	vtkSmartPointer<vtkDataSetMapper> lineMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	lineMapper->SetInputConnection(line->GetOutputPort());
	vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
	lineActor->SetMapper(lineMapper);
	lineActor->GetProperty()->SetColor(1, 0, 0);
	lineActor->GetProperty()->SetLineWidth(3);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(lineActor);
	renderer->SetBackground(1, 1, 1);
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);
	renWin->SetSize(300, 300);
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(style);


	vtkSmartPointer<vtkAnimationScene> scene = vtkAnimationScene::New();
	scene->SetModeToSequence();
	scene->SetFrameRate(10);
	scene->SetStartTime(0);
	scene->SetEndTime(100);


	vtkCustomAnimationCue* cue1 = vtkCustomAnimationCue::New();
	cue1->line = line;
	cue1->points = points;
	cue1->RenWin = renWin;  //对于类中需要操作的对象进行赋值
	cue1->SetTimeModeToNormalized();
	cue1->SetStartTime(0);
	cue1->SetEndTime(1.0);
	scene->AddCue(cue1);
	scene->Play();
	scene->Stop();
	cue1->Delete();



	iren->Initialize();
	iren->Start();
	return 0;

}
*/


//动画显示，连续显示文件测试，未完成，（for循环导致线程堵塞）
/*
#include <vtkAppendFilter.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <algorithm>
#include <array>
#include <string>
#include <vtkScalarBarActor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkStdString.h>


namespace {
	vtkSmartPointer<vtkUnstructuredGrid>
		ReadUnstructuredGrid(std::string const& fileName);
}

int main(int argc, char* argv[])
{
	// Vis Pipeline
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(640, 480);
	renderWindow->AddRenderer(renderer);
	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);

	renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
	renderer->UseHiddenLineRemovalOn();

	renderer->GetActiveCamera()->Azimuth(45);
	renderer->GetActiveCamera()->Elevation(45);

	//坐标轴
	//vtkNew<vtkAxesActor> axes;
	//vtkNew<vtkOrientationMarkerWidget> widget;
	//double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
	//colors->GetColor("Carrot", rgba);
	//widget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
	//widget->SetOrientationMarker(axes);
	//widget->SetInteractor(interactor);
	//widget->SetViewport(0.0, 0.0, 0.4, 0.4);
	//widget->SetEnabled(1);
	//widget->InteractiveOn();

	vtkNew<vtkUnstructuredGridReader> reader;
	vtkNew<vtkLookupTable> lut1;
	vtkNew<vtkDataSetMapper> mapper;
	vtkNew<vtkScalarBarActor> scalarbar;
	vtkNew<vtkActor> actor;


	//循环读入多个.vtk文件
	vtkSmartPointer<vtkStringArray > fileArray =
		vtkSmartPointer<vtkStringArray >::New();
	char fileName[128];
	for (int i = 0; i < 11; i++) //几个图像就循环几次
	{
		sprintf_s(fileName, "Test%02d.vtk", i+1);
		vtkStdString fileStr(fileName);
		fileArray->InsertNextValue(fileStr);
		//fileArray->SetValue(i, fileStr);
	}
	for (int j = 0; j < 10; ++j)
	{
		for (int i = 0; i < 11; i++)
		{
			std::string name;
			name = fileArray->GetValue(i);
			//vtkStdString setFileName = fileArray;

			reader->SetFileName(name.c_str());
			reader->Update();

			//std::cout << "Loading: " << argv[1] << std::endl;
			//auto unstructuredGrid = ReadUnstructuredGrid(std::string(argv[1]));
			auto unstructuredGrid = reader->GetOutputPort();
			auto unstructuredgrid = reader->GetOutput();


			lut1->SetHueRange(0.5, 0.833);
			// Visualize

			//mapper->SetInputData(unstructuredGrid);
			mapper->SetInputConnection(unstructuredGrid);
			//mapper->ScalarVisibilityOff();
			mapper->SetScalarRange(unstructuredgrid->GetScalarRange());
			mapper->SetLookupTable(lut1);
			mapper->SetColorModeToMapScalars();


			scalarbar->SetLookupTable(mapper->GetLookupTable());
			//scalarbar->SetTitle(curvaturesfilter->GetOutput()->GetPointData()->GetScalars()->GetName());
			scalarbar->SetNumberOfLabels(5);
			renderer->AddActor2D(scalarbar);





			actor->SetMapper(mapper);
			//actor->SetBackfaceProperty(backProp);
			//actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
			//actor->GetProperty()->SetSpecular(.3);
			//actor->GetProperty()->SetSpecularPower(30);
			//actor->GetProperty()->EdgeVisibilityOn();
			renderer->AddActor(actor);
			renderer->ResetCamera();
			renderer->Modified();
			//renderWindow->Render();
			interactor->Initialize();
			//interactor->Start();
			interactor->Render();//renderWindow以及调用render()，所以可以不要
			//Delay(0.5 * 1000);   //延时0.5秒

		}
	}
	interactor->Start();
	return EXIT_SUCCESS;
}
*/

/*
#include <vtkActor.h>
#include <vtkAffineRepresentation2D.h>
#include <vtkAffineWidget.h>
#include <vtkAppendPolyData.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>

namespace {
	class vtkAffineCallback : public vtkCommand
	{
	public:
		static vtkAffineCallback* New()
		{
			return new vtkAffineCallback;
		}
		virtual void Execute(vtkObject* caller, unsigned long, void*);
		vtkAffineCallback() : Actor(0), AffineRep(0)
		{
			this->Transform = vtkTransform::New();
		}
		~vtkAffineCallback()
		{
			this->Transform->Delete();
		}
		vtkActor* Actor;
		vtkAffineRepresentation2D* AffineRep;
		vtkTransform* Transform;
	};
} // namespace

int main(int, char* [])
{
	vtkNew<vtkNamedColors> colors;

	// Create two spheres: a larger one and a smaller one on top of the larger one
	// to show a reference point while rotating
	vtkNew<vtkSphereSource> sphereSource;
	sphereSource->Update();

	vtkNew<vtkSphereSource> sphereSource2;
	sphereSource2->SetRadius(0.075);
	sphereSource2->SetCenter(0, 0.5, 0);
	sphereSource2->Update();

	// Append the two spheres into one vtkPolyData
	vtkNew<vtkAppendPolyData> append;
	append->AddInputConnection(sphereSource->GetOutputPort());
	append->AddInputConnection(sphereSource2->GetOutputPort());

	// Create a plane centered over the larger sphere with 4x4 sub sections
	vtkNew<vtkPlaneSource> planeSource;
	planeSource->SetXResolution(4);
	planeSource->SetYResolution(4);
	planeSource->SetOrigin(-1, -1, 0);
	planeSource->SetPoint1(1, -1, 0);
	planeSource->SetPoint2(-1, 1, 0);

	// Create a mapper and actor for the plane: show it as a wireframe
	vtkNew<vtkPolyDataMapper> planeMapper;
	planeMapper->SetInputConnection(planeSource->GetOutputPort());
	vtkNew<vtkActor> planeActor;
	planeActor->SetMapper(planeMapper);
	planeActor->GetProperty()->SetRepresentationToWireframe();
	planeActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

	// Create a mapper and actor for the spheres
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(append->GetOutputPort());
	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);

	// Create a renderer and render window
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);
	renderWindow->SetWindowName("AffineWidget");

	renderer->AddActor(actor);
	renderer->AddActor(planeActor);
	renderer->GradientBackgroundOn();
	renderer->SetBackground(colors->GetColor3d("LightSkyBlue").GetData());
	renderer->SetBackground2(colors->GetColor3d("MidnightBlue").GetData());

	// Create an interactor
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);
	dynamic_cast<vtkInteractorStyleSwitch*>(
		renderWindowInteractor->GetInteractorStyle())
		->SetCurrentStyleToTrackballCamera();

	// Create an affine widget to manipulate the actor
	// the widget currently only has a 2D representation and therefore applies
	// transforms in the X-Y plane only
	vtkNew<vtkAffineWidget> affineWidget;
	affineWidget->SetInteractor(renderWindowInteractor);
	affineWidget->CreateDefaultRepresentation();
	dynamic_cast<vtkAffineRepresentation2D*>(affineWidget->GetRepresentation())
		->PlaceWidget(actor->GetBounds());

	vtkNew<vtkAffineCallback> affineCallback;
	affineCallback->Actor = actor;
	affineCallback->AffineRep = dynamic_cast<vtkAffineRepresentation2D*>(
		affineWidget->GetRepresentation());

	affineWidget->AddObserver(vtkCommand::InteractionEvent, affineCallback);
	affineWidget->AddObserver(vtkCommand::EndInteractionEvent, affineCallback);

	renderWindow->Render();
	renderWindowInteractor->Initialize();
	renderWindow->Render();
	affineWidget->On();

	// begin mouse interaction
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}

namespace {
	void vtkAffineCallback::Execute(vtkObject*, unsigned long vtkNotUsed(event),
		void*)
	{
		this->AffineRep->GetTransform(this->Transform);
		this->Actor->SetUserTransform(this->Transform);
	}
} // namespace
*/



//球移动动画，vtkcommond的时间命令类扩展，鼠标不被阻塞，
/*
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>

namespace {
	class vtkTimerCallback2 : public vtkCommand
	{
	public:
		vtkTimerCallback2() = default;

		static vtkTimerCallback2* New()
		{
			vtkTimerCallback2* cb = new vtkTimerCallback2;
			cb->TimerCount = 0;
			return cb;
		}

		virtual void Execute(vtkObject* caller, unsigned long eventId,
			void* vtkNotUsed(callData))
		{
			if (vtkCommand::TimerEvent == eventId)
			{
				++this->TimerCount;
			}
			std::cout << this->TimerCount << std::endl;
			actor->SetPosition(this->TimerCount, this->TimerCount, 0);
			if (this->TimerCount < this->maxCount)
			{

				vtkRenderWindowInteractor* iren =
					dynamic_cast<vtkRenderWindowInteractor*>(caller);
				iren->GetRenderWindow()->Render();
			}
			else
			{
				vtkRenderWindowInteractor* iren =
					dynamic_cast<vtkRenderWindowInteractor*>(caller);
				if (this->timerId > -1)
				{
					iren->DestroyTimer(this->timerId);
				}
			}
		}

	private:
		int TimerCount = 0;

	public:
		vtkActor* actor = nullptr;
		int timerId = 0;
		int maxCount = -1;
	};
} // namespace

int main(int, char* [])
{
	vtkNew<vtkNamedColors> colors;

	// Create a sphere
	vtkNew<vtkSphereSource> sphereSource;
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(2.0);
	sphereSource->SetPhiResolution(30);
	sphereSource->SetThetaResolution(30);

	// Create a mapper and actor
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(sphereSource->GetOutputPort());
	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetSpecular(0.6);
	actor->GetProperty()->SetSpecularPower(30);
	actor->GetProperty()->SetColor(colors->GetColor3d("Peacock").GetData());

	// Create a renderer, render window, and interactor
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);
	renderWindow->SetWindowName("Animation");

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actor to the scene
	renderer->AddActor(actor);
	renderer->SetBackground(colors->GetColor3d("MistyRose").GetData());

	// Render and interact
	renderWindow->Render();
	renderer->GetActiveCamera()->Zoom(0.8);
	renderWindow->Render();

	// Initialize must be called prior to creating timer events.
	renderWindowInteractor->Initialize();

	// Sign up to receive TimerEvent
	vtkNew<vtkTimerCallback2> cb;
	cb->actor = actor;
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, cb);

	int timerId = renderWindowInteractor->CreateRepeatingTimer(500);
	std::cout << "timerId: " << timerId << std::endl;
	// Destroy the timer when maxCount is reached.
	cb->maxCount = 9;
	cb->timerId = timerId;
	// Start the interaction and timer
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
*/

//动画显示，连续显示文件测试，commond Time事件

//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkLookupTable.h>
//#include <algorithm>
//#include <array>
//#include <string>
//#include <vtkScalarBarActor.h>
//#include <vtkAxesActor.h>
//#include <vtkOrientationMarkerWidget.h>
//#include <vtkStdString.h>
//
//static unsigned int fileNumber = 0;
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//namespace {
//	class vtkTimerCallback2 : public vtkCommand
//	{
//	public:
//		vtkTimerCallback2() = default;
//
//		static vtkTimerCallback2* New()
//		{
//			vtkTimerCallback2* cb = new vtkTimerCallback2;
//			//fileNumber = 0;
//			return cb;
//		}
//
//		virtual void Execute(vtkObject* caller, unsigned long eventId,
//			void* vtkNotUsed(callData))
//		{
//			if (vtkCommand::TimerEvent == eventId)
//			{
//				//++this->TimerCount;
//				++fileNumber;
//			}
//
//			//actor->SetPosition(this->TimerCount, this->TimerCount, 0);
//				std::string name;
//				name = array->GetValue(fileNumber);
//				//vtkStdString setFileName = fileArray;
//
//				reader->SetFileName(name.c_str());
//				reader->Update();
//
//				auto unstructuredGrid = reader->GetOutputPort();
//				auto unstructuredgrid = reader->GetOutput();
//
//
//				//lut1->SetHueRange(0.5, 0.833);
//
//
//				//mapper->SetInputData(unstructuredGrid);
//				mapper->SetInputConnection(unstructuredGrid);
//				//mapper->ScalarVisibilityOff();
//				mapper->SetScalarRange(unstructuredgrid->GetScalarRange());
//				mapper->SetLookupTable(lut);
//				mapper->SetColorModeToMapScalars();
//				actor->SetMapper(mapper);
//
//			//if (this->TimerCount < this->maxCount)
//			if (fileNumber < 11)
//			{
//				
//
//				vtkRenderWindowInteractor* iren =
//					dynamic_cast<vtkRenderWindowInteractor*>(caller);
//				iren->GetRenderWindow()->Render();
//			}
//			else
//			{
//				fileNumber = 0;
//				vtkRenderWindowInteractor* iren =
//					dynamic_cast<vtkRenderWindowInteractor*>(caller);
//				iren->GetRenderWindow()->Render();
//				/*if (this->timerId > -1)
//				{
//					//iren->DestroyTimer(this->timerId);
//					fileNumber = 0;
//				}*/
//			}
//		}
//
//	private:
//		int TimerCount = 0;
//
//	public:
//		vtkActor* actor = nullptr;
//		int timerId = 0;
//		int maxCount = -1;
//		vtkStringArray* array = nullptr;
//		vtkUnstructuredGridReader* reader = nullptr;
//		vtkMapper* mapper = nullptr;
//		vtkLookupTable* lut;
//	};
//} // namespace
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//
//	/*
//	坐标轴
//	vtkNew<vtkAxesActor> axes;
//	vtkNew<vtkOrientationMarkerWidget> widget;
//	double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
//	colors->GetColor("Carrot", rgba);
//	widget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
//	widget->SetOrientationMarker(axes);
//	widget->SetInteractor(interactor);
//	widget->SetViewport(0.0, 0.0, 0.4, 0.4);
//	widget->SetEnabled(1);
//	widget->InteractiveOn();
//	*/
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	vtkNew<vtkLookupTable> lut1;
//	vtkNew<vtkDataSetMapper> mapper;
//	vtkNew<vtkScalarBarActor> scalarbar;
//	vtkNew<vtkActor> actor;
//
//
//	//循环读入多个.vtk文件
//	vtkSmartPointer<vtkStringArray > fileArray =
//		vtkSmartPointer<vtkStringArray >::New();
//	char fileName[128];
//	for (int i = 0; i < 11; i++) //几个图像就循环几次
//	{
//		sprintf_s(fileName, "Test%02d.vtk", i + 1);
//		vtkStdString fileStr(fileName);
//		fileArray->InsertNextValue(fileStr);
//		//fileArray->SetValue(i, fileStr);
//	}
//
//
//	std::string name;
//	name = fileArray->GetValue(fileNumber);
//	//vtkStdString setFileName = fileArray;
//
//	reader->SetFileName(name.c_str());
//	reader->Update();
//
//	//std::cout << "Loading: " << argv[1] << std::endl;
//	//auto unstructuredGrid = ReadUnstructuredGrid(std::string(argv[1]));
//	
//	auto unstructuredGrid = reader->GetOutputPort();
//	auto unstructuredgrid = reader->GetOutput();
//	
//	
//	//lut1->SetHueRange(0.5, 0.833);
//	
//	
//	//mapper->SetInputData(unstructuredGrid);
//	mapper->SetInputConnection(unstructuredGrid);
//	//mapper->ScalarVisibilityOff();
//	mapper->SetScalarRange(unstructuredgrid->GetScalarRange());
//	mapper->SetLookupTable(lut1);
//	mapper->SetColorModeToMapScalars();
//
//
//	
//	/*scalarbar->SetLookupTable(mapper->GetLookupTable());
//	//scalarbar->SetTitle(curvaturesfilter->GetOutput()->GetPointData()->GetScalars()->GetName());
//	scalarbar->SetNumberOfLabels(5);
//	renderer->AddActor2D(scalarbar);
//	*/
//
//
//
//
//	actor->SetMapper(mapper);
//	//actor->SetBackfaceProperty(backProp);
//	//actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
//	//actor->GetProperty()->SetSpecular(.3);
//	//actor->GetProperty()->SetSpecularPower(30);
//	//actor->GetProperty()->EdgeVisibilityOn();
//	renderer->AddActor(actor);
//	renderer->ResetCamera();
//	//renderer->Modified();
//	renderWindow->Render();
//	//renderWindow->Render();
//	interactor->Initialize();
//
//	//interactor->Start();
//	//interactor->Render();//renderWindow以及调用render()，所以可以不要
//	//Delay(0.5 * 1000);   //延时0.5秒
//
//
//
//	vtkNew<vtkTimerCallback2> cb;
//	cb->actor = actor;
//	interactor->AddObserver(vtkCommand::TimerEvent, cb);
//
//	int timerId = interactor->CreateRepeatingTimer(500);
//	std::cout << "timerId: " << timerId << std::endl;
//
//	// Destroy the timer when maxCount is reached.
//	cb->maxCount = 9;
//	cb->timerId = timerId;
//	cb->array = fileArray;
//	cb->reader = reader;
//	cb->mapper = mapper;
//	cb->lut = lut1;
//	// Start the interaction and timer
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}


//左键选择actor，高亮显示

//#include <vtkVersion.h>
//#include <vtkSmartPointer.h>
//#include <vtkMath.h>
//#include <vtkActor.h>
//#include <vtkProperty.h>
//#include <vtkSphereSource.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkObjectFactory.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPropPicker.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//
//// Handle mouse events
//class MouseInteractorHighLightActor : public vtkInteractorStyleTrackballCamera
//{
//public:
//	static MouseInteractorHighLightActor* New();
//	vtkTypeMacro(MouseInteractorHighLightActor, vtkInteractorStyleTrackballCamera);
//
//	MouseInteractorHighLightActor()
//	{
//		LastPickedActor = NULL;
//		LastPickedProperty = vtkProperty::New();
//	}
//	virtual ~MouseInteractorHighLightActor()
//	{
//		LastPickedProperty->Delete();
//	}
//	virtual void OnLeftButtonDown()
//	{
//		int* clickPos = this->GetInteractor()->GetEventPosition();
//
//		// Pick from this location.
//		vtkSmartPointer<vtkPropPicker>  picker =
//			vtkSmartPointer<vtkPropPicker>::New();
//		picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
//
//		// If we picked something before, reset its property
//		if (this->LastPickedActor)
//		{
//			this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
//		}
//		this->LastPickedActor = picker->GetActor();
//		if (this->LastPickedActor)
//		{
//			// Save the property of the picked actor so that we can
//			// restore it next time
//			this->LastPickedProperty->DeepCopy(this->LastPickedActor->GetProperty());
//			// Highlight the picked actor by changing its properties
//			this->LastPickedActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
//			this->LastPickedActor->GetProperty()->SetDiffuse(1.0);
//			this->LastPickedActor->GetProperty()->SetSpecular(0.0);
//		}
//
//		// Forward events
//		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
//	}
//
//private:
//	vtkActor* LastPickedActor;
//	vtkProperty* LastPickedProperty;
//};
//
//vtkStandardNewMacro(MouseInteractorHighLightActor);
//
//// Execute application.
//int main(int argc, char* argv[])
//{
//	int numberOfSpheres = 10;
//	if (argc > 1)
//	{
//		numberOfSpheres = atoi(argv[1]);
//	}
//	// A renderer and render window
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	// An interactor
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Set the custom type to use for interaction.
//	vtkSmartPointer<MouseInteractorHighLightActor> style =
//		vtkSmartPointer<MouseInteractorHighLightActor>::New();
//	style->SetDefaultRenderer(renderer);
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	for (int i = 0; i < numberOfSpheres; ++i)
//	{
//		vtkSmartPointer<vtkSphereSource> source =
//			vtkSmartPointer<vtkSphereSource>::New();
//		double x, y, z, radius;
//		x = vtkMath::Random(-5, 5);
//		y = vtkMath::Random(-5, 5);
//		z = vtkMath::Random(-5, 5);
//		radius = vtkMath::Random(.5, 1.0);
//		source->SetRadius(radius);
//		source->SetCenter(x, y, z);
//		source->SetPhiResolution(11);
//		source->SetThetaResolution(21);
//		vtkSmartPointer<vtkPolyDataMapper> mapper =
//			vtkSmartPointer<vtkPolyDataMapper>::New();
//		mapper->SetInputConnection(source->GetOutputPort());
//		vtkSmartPointer<vtkActor> actor =
//			vtkSmartPointer<vtkActor>::New();
//		actor->SetMapper(mapper);
//		double r, g, b;
//		r = vtkMath::Random(.4, 1.0);
//		g = vtkMath::Random(.4, 1.0);
//		b = vtkMath::Random(.4, 1.0);
//		actor->GetProperty()->SetDiffuseColor(r, g, b);
//		actor->GetProperty()->SetDiffuse(.8);
//		actor->GetProperty()->SetSpecular(.5);
//		actor->GetProperty()->SetSpecularColor(1.0, 1.0, 1.0);
//		actor->GetProperty()->SetSpecularPower(30.0);
//		renderer->AddActor(actor);
//	}
//
//	renderer->SetBackground(.3, .4, .5);
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Initialize();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//摆线，不能交互

//#include "vtkSmartPointer.h"
//#include "vtkRenderWindow.h"
//#include "vtkPoints.h"
//#include "vtkParametricSpline.h"
//#include "vtkSetGet.h"
//#include "vtkObjectFactory.h"
//#include "vtkParametricFunctionSource.h"
//#include "vtkPolyDataMapper.h"
//#include "vtkActor.h"
//#include "vtkRenderer.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkInteractorStyleTrackballCamera.h"
//#include "vtkAnimationScene.h"
//#include "vtkAnimationCue.h"
//#include "vtkProperty.h"
//
//class Spline_Cue : public vtkAnimationCue
//{
//public:
//	vtkTypeMacro(Spline_Cue, vtkAnimationCue);
//	static Spline_Cue* New();
//
//	vtkSmartPointer<vtkRenderWindow> RenWin;
//	vtkSmartPointer<vtkPoints> points;
//	vtkSmartPointer<vtkParametricSpline> spline;
//protected:
//	Spline_Cue()
//	{
//		this->RenWin = nullptr;
//		this->points = nullptr;
//		this->spline = nullptr;
//	}
//	virtual void TickInternal(double currenttime, double deltatime, double clocktime)
//	{//在每一帧调用的重新绘制函数
//		double x = 10.0, z = 0;
//		double y = 4 * sin(2 * 3.14 * currenttime);
//		this->points->SetPoint(2, x, y, z);
//		this->spline->Modified();
//		this->RenWin->Render();//重新渲染
//	}
//};
//vtkStandardNewMacro(Spline_Cue);
//
//int main()
//{
//	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//	points->InsertPoint(0, 0.0, 0.0, 0.0);
//	points->InsertPoint(1, 5.0, 0, 0.0);
//	points->InsertPoint(2, 10.0, 10, 0.0);
//
//	vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
//	spline->SetPoints(points);
//	spline->ClosedOff();//不需要闭合
//
//	vtkSmartPointer<vtkParametricFunctionSource> splineSource =
//		vtkSmartPointer<vtkParametricFunctionSource>::New();
//	splineSource->SetParametricFunction(spline);
//
//	vtkSmartPointer<vtkPolyDataMapper> splineMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	splineMapper->SetInputConnection(splineSource->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();//创建角色
//	lineActor->SetMapper(splineMapper);
//	lineActor->GetProperty()->SetColor(1, 0, 0);
//	lineActor->GetProperty()->SetLineWidth(3);
//
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();//创建渲染场景
//	renderer->AddActor(lineActor);
//	renderer->SetBackground(0, 0, 0);
//	vtkSmartPointer<vtkRenderWindow> renWin =
//		vtkSmartPointer<vtkRenderWindow>::New();//创建显示窗口
//	renWin->AddRenderer(renderer);
//	renWin->SetSize(600, 600);
//	vtkSmartPointer<vtkRenderWindowInteractor> iren =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();//创建交互器
//	iren->SetRenderWindow(renWin);
//	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
//		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
//	iren->SetInteractorStyle(style);
//
//	vtkSmartPointer<vtkAnimationScene> scene = vtkSmartPointer<vtkAnimationScene>::New();//创建动画场景
//	scene->SetModeToSequence();//设置顺序播放模式
////	scene->SetModeToRealTime();//设置实时播放模式
//	scene->SetFrameRate(30);//设置帧率,单位时间内渲染的帧数
//	scene->SetStartTime(0);//动画开始时间
//	scene->SetEndTime(10);//动画结束时间
//
//	vtkSmartPointer<Spline_Cue> cue1 = vtkSmartPointer<Spline_Cue>::New();//创建动画实例
//	cue1->RenWin = renWin;//对于类中需要操作的对象进行赋值  
//	cue1->points = points;
//	cue1->spline = spline;
//
//	cue1->SetTimeModeToNormalized();//按照场景时间标准化实例的动画时间：0对应动画场景的开始，1对应结束
//	cue1->SetStartTime(0);
//	cue1->SetEndTime(1.0);
//	scene->AddCue(cue1);//添加动画实例到动画场景中
//	scene->SetLoop(1);//设置循环播放模式
//	scene->Play();//动画场景开始播放
//	scene->Stop();//停止播放
//
//	iren->Initialize();
//	iren->Start();
//	return 0;
//}

//写入点文件

//#include<vtkSmartPointer.h>
//#include<vtkPoints.h>
//#include<vtkWriter.h>
//#include<vtkPolyDataWriter.h>
//#include<vtkPolyData.h>
//#include<vtkCellArray.h>
//
//int main(int argc, char* argv[])
//{
	//vtkSmartPointer<vtkPoints>points =
	//	vtkSmartPointer<vtkPoints>::New();
	//points->InsertNextPoint(1.0, 0.0, 0.0);
	//points->InsertNextPoint(0.0, 0.0, 0.0);
	//points->InsertNextPoint(0.0, 1.0, 0.0);

//vtkSmartPointer<vtkPolyData>polydata =
//	vtkSmartPointer<vtkPolyData>::New();
//polydata->SetPoints(points);
//
//vtkSmartPointer<vtkPolyDataWriter>writer =
//	vtkSmartPointer<vtkPolyDataWriter>::New();
//writer->SetFileName("triangle.vtk");
//writer->SetInputData(polydata);
//writer->Write();
//return 0;
//
//	double X[3] = { 1.0,0.0,0.0 };
//	double Y[3] = { 0.0,0.0,1.0 };
//	double Z[3] = {	0.0,0.0,0.0 };
//	
//	vtkSmartPointer<vtkPoints>points =
//		vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkCellArray>vertices =
//		vtkSmartPointer<vtkCellArray>::New();
//
//	for (unsigned int i = 0; i < 3; ++i)
//	{
//		vtkIdType pid[1];
//		pid[0] = points->InsertNextPoint(X[i], Y[i], Z[i]);
//		vertices->InsertNextCell(1, pid);
//	}
//
//	vtkSmartPointer<vtkPolyData>polydata =
//		vtkSmartPointer<vtkPolyData>::New();
//
//	polydata->SetPoints(points);
//	polydata->SetVerts(vertices);
//
//	vtkSmartPointer<vtkPolyDataWriter>writer =
//		vtkSmartPointer<vtkPolyDataWriter>::New();
//	writer->SetFileName("TriangleVerts.vtk");
//	writer->SetInputData(polydata);
//	writer->Write();
//	return 0;
//}

//随机高度的平面

//#include<vtkSmartPointer.h>
//#include<vtkPoints.h>
//#include<vtkPolyData.h>
//#include<vtkDelaunay2D.h>
//#include<vtkDataSetMapper.h>
//#include<vtkActor.h>
//#include<vtkRenderer.h>
//#include<vtkRenderWindow.h>
//#include<vtkRenderWindowInteractor.h>
//
//
//int main(int argc, char* argv[])
//{
//	vtkSmartPointer<vtkPoints> points =
//		vtkSmartPointer<vtkPoints>::New();
//	unsigned int gridSize = 10;
//	for (unsigned int x = 0; x < gridSize; ++x)
//	{
//		for (unsigned int y = 0; y < gridSize; ++y)
//		{
//			points->InsertNextPoint(x, y, vtkMath::Random(0.0, 3.0));
//		}
//	}
//	vtkSmartPointer<vtkPolyData>polydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	polydata->SetPoints(points);
//
//	vtkSmartPointer<vtkDelaunay2D>delaunay =
//		vtkSmartPointer<vtkDelaunay2D>::New();
//	delaunay->SetInputData(polydata);
//	delaunay->Update();
//
//	vtkSmartPointer<vtkDataSetMapper>mapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	mapper->SetInputConnection(delaunay->GetOutputPort());
//
//	vtkSmartPointer<vtkActor>actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	vtkSmartPointer<vtkRenderer>renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	renderer->AddActor(actor);
//	renderer->SetBackground(1.0, 0.3, 0.2);
//	renderer->SetViewport(0.0, 0.0, 0.5, 0.5);
//
//	vtkSmartPointer<vtkRenderWindow>renWin =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renWin->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor>iren =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	iren->Initialize();
//	iren->SetRenderWindow(renWin);
//	iren->Start();
//	renWin->Render();
//}


//剖切

//#include <vtkSmartPointer.h>
//#include <vtkClipDataSet.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkPlane.h>
//#include <vtkTransform.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkCellTypes.h>
//#include <vtkDataSetMapper.h>
//#include <vtkLookupTable.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkRenderer.h>
//#include <vtkNamedColors.h>

//int main(int argc, char* argv[])
//{
//
//	// Create the reader for the data.
//
//
//	//auto reader =
//	//	vtkSmartPointer<vtkUnstructuredGridReader>::New();
//	auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
//	reader->SetFileName("vtuTest.vtu");
//	reader->Update();
//
//	double bounds[6];
//	reader->GetOutput()->GetBounds(bounds);
//	double center[3];
//	reader->GetOutput()->GetCenter(center);
//
//	auto colors =
//		vtkSmartPointer<vtkNamedColors>::New();
//	auto renderer = vtkSmartPointer<vtkRenderer>::New();
//	renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	auto renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetSize(640, 480);
//
//	auto interactor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	interactor->SetRenderWindow(renderWindow);
//
//	double xnorm[3] = { 1.0, 0, 0 };//切割面法线方向
//
//	//创建切割平面，设置中点和法线方向
//	auto clipPlane = vtkSmartPointer<vtkPlane>::New();
//	clipPlane->SetOrigin(reader->GetOutput()->GetCenter());
//	//clipPlane->SetOrigin(0.5,0.5,1);
//	clipPlane->SetNormal(xnorm);
//
//	//创建切割面数据
//	auto clipper =
//		vtkSmartPointer<vtkClipDataSet>::New();
//	clipper->SetClipFunction(clipPlane);
//	clipper->SetInputData(reader->GetOutput());
//	clipper->SetValue(0.0);
//	clipper->GenerateClippedOutputOn();
//	clipper->Update();
//
//	//创建切割面映射
//	auto insideMapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	insideMapper->SetInputData(clipper->GetOutput());
//	//insideMapper->ScalarVisibilityOff();
//	
//
//	vtkNew<vtkLookupTable> lut1;
//	lut1->SetHueRange(0.5, 0.833);
//	// Visualize
//
//	//mapper->ScalarVisibilityOff();
//	insideMapper->SetScalarRange(reader->GetOutput()->GetScalarRange());
//	insideMapper->SetLookupTable(lut1);
//	insideMapper->SetColorModeToMapScalars();
//
//	auto insideActor =
//		vtkSmartPointer<vtkActor>::New();
//	insideActor->SetMapper(insideMapper);
//	//insideActor->GetProperty()->SetDiffuseColor(colors->GetColor3d("banana").GetData());
//	insideActor->GetProperty()->SetAmbient(.3);
//	//设置是否显示网格
//	insideActor->GetProperty()->EdgeVisibilityOn();
//
//
//
//
//
//	/*auto clippedMapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	clippedMapper->SetInputData(clipper->GetClippedOutput());
//	clippedMapper->ScalarVisibilityOff();
//
//	auto clippedActor =
//		vtkSmartPointer<vtkActor>::New();
//	clippedActor->SetMapper(clippedMapper);
//	clippedActor->GetProperty()->SetDiffuseColor(colors->GetColor3d("tomato").GetData());
//	insideActor->GetProperty()->SetAmbient(.3);
//	clippedActor->GetProperty()->EdgeVisibilityOn();*/
//
//	// Create transforms to make a better visualization
//	auto insideTransform = vtkSmartPointer<vtkTransform>::New();
//	insideTransform->Translate(-(bounds[1] - bounds[0]) * .75, 0, 0);
//	insideTransform->Translate(center[0], center[1], center[2]);
//	insideTransform->RotateY(-120.0);
//	insideTransform->Translate(-center[0], -center[1], -center[2]);
//	insideActor->SetUserTransform(insideTransform);
//
//	//auto clippedTransform = vtkSmartPointer<vtkTransform>::New();
//	//clippedTransform->Translate((bounds[1] - bounds[0]) * .75, 0, 0);
//	//clippedTransform->Translate(center[0], center[1], center[2]);
//	//clippedTransform->RotateY(60.0);
//	//clippedTransform->Translate(-center[0], -center[1], -center[2]);
//	//clippedActor->SetUserTransform(clippedTransform);
//
//	//renderer->AddViewProp(clippedActor);
//	renderer->AddViewProp(insideActor);
//
//	renderer->ResetCamera();
//	renderer->GetActiveCamera()->Dolly(1.4);
//	renderer->ResetCameraClippingRange();
//	renderWindow->Render();
//	renderWindow->SetWindowName("ClipUnstructuredGridWithPlane2");
//	renderWindow->Render();
//
//	interactor->Start();
//
//	
//	return EXIT_SUCCESS;
//}

//vtkPlaneWidget网格面

//#include <vtkPlaneWidget.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//
//int main(int, char* [])
//{
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	vtkSmartPointer<vtkPlaneWidget> planeWidget =
//		vtkSmartPointer<vtkPlaneWidget>::New();
//	planeWidget->SetInteractor(renderWindowInteractor);
//
//	planeWidget->On();
//
//	renderWindowInteractor->Initialize();
//
//	renderer->ResetCamera();
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//剖切，切块。改为了提取表面一层进行渲染，缩短渲染时间

//#include <vtkSmartPointer.h>
//
//#include <vtkXMLPolyDataReader.h>
//#include <vtkSphereSource.h>
//#include <vtkClipPolyData.h>
//#include <vtkPlane.h>
//
//#include <vtkCommand.h>
//#include <vtkImplicitPlaneWidget2.h>
//#include <vtkImplicitPlaneRepresentation.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNew.h>
//#include <vtkLookupTable.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkClipDataSet.h>
//#include <vtkCutter.h>
//#include <vtkMapper.h>
//
//// Callback for the interaction
//// This does the actual work: updates the vtkPlane implicit function.
//// This in turn causes the pipeline to update and clip the object.
//class vtkIPWCallback : public vtkCommand
//{
//public:
//	static vtkIPWCallback* New()
//	{
//		return new vtkIPWCallback;
//	}
//	virtual void Execute(vtkObject* caller, unsigned long, void*)
//	{
//		vtkImplicitPlaneWidget2* planeWidget =
//			reinterpret_cast<vtkImplicitPlaneWidget2*>(caller);
//		vtkImplicitPlaneRepresentation* rep =
//			reinterpret_cast<vtkImplicitPlaneRepresentation*>(planeWidget->GetRepresentation());
//		rep->GetPlane(this->Plane);
//	}
//	vtkIPWCallback() :Plane(0), Actor(0) {}
//	vtkPlane* Plane;
//	vtkActor* Actor;
//
//};
//
//int main(int argc, char* argv[])
//{
//
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	vtkSmartPointer<vtkUnstructuredGridReader> reader =
//		vtkSmartPointer<vtkUnstructuredGridReader>::New();
//	//vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
//	//	vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
//	reader->SetFileName("1.vtk");
//	reader->Update();
//
//	auto unstructuredgrid = reader->GetOutput();
//
//	vtkSmartPointer<vtkDataSetMapper>mapper0 =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	mapper0->SetInputData(unstructuredgrid);
//	//mapper->SetScalarRange(this->getScalarRange(scalarIte));
//
//	//零件actor
//	vtkSmartPointer<vtkActor>opacityActor =
//		vtkSmartPointer<vtkActor>::New();
//	opacityActor->SetMapper(mapper0);
//	opacityActor->GetProperty()->SetOpacity(0.1);
//	// 添加renderer的元素
//	renderer->AddActor(opacityActor);
//
//	vtkNew<vtkLookupTable> lut1;
//	lut1->SetHueRange(0.4, 0.9);
//
//	//vtkNew<vtkMapper>sourceMapper;
//	//vtkSmartPointer<vtkMapper>sourceMapper =
//	//	vtkSmartPointer<vtkMapper>::New();
//	vtkMapper* sourceMapper = renderer->GetActors()->GetLastActor()->GetMapper();
//
//	//抽出表面
//	vtkSmartPointer<vtkDataSetSurfaceFilter> ds_filter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//	ds_filter->SetInputData(sourceMapper->GetInput());
//	ds_filter->Update();
//
//	//polydata数据
//	vtkSmartPointer<vtkPolyData> polyData = ds_filter->GetOutput();
//	//剖切平面
//	vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
//	plane->SetOrigin(reader->GetOutput()->GetCenter());
//
//	//clipper,出去平面剩下的部分，很显然由于ds_filter把表面提取出来了，因此这个只是一个无盖的壳
//	vtkSmartPointer<vtkClipPolyData> clipperData = vtkSmartPointer<vtkClipPolyData>::New();
//	clipperData->SetInputData(polyData);
//	clipperData->SetClipFunction(plane);
//	clipperData->Modified();
//	//cutter,只有一个平面
//	vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
//	cutter->SetInputData(sourceMapper->GetInput());
//	cutter->SetCutFunction(plane);
//	cutter->Modified();
//	//clip mapper
//	vtkSmartPointer<vtkPolyDataMapper> clip_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	clip_mapper->SetLookupTable(lut1);
//	clip_mapper->SetScalarRange(polyData->GetScalarRange());
//	clip_mapper->SetInputConnection(clipperData->GetOutputPort());
//	//cut mapper
//	vtkSmartPointer<vtkPolyDataMapper> cut_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	cut_mapper->SetLookupTable(lut1);
//	cut_mapper->SetScalarRange(sourceMapper->GetInput()->GetScalarRange());
//	cut_mapper->SetInputConnection(cutter->GetOutputPort());
//	//actor
//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(clip_mapper);
//	vtkSmartPointer<vtkActor> cut_actor = vtkSmartPointer<vtkActor>::New();
//	cut_actor->SetMapper(cut_mapper);
//
//	vtkSmartPointer<vtkIPWCallback> myCallback =
//		vtkSmartPointer<vtkIPWCallback>::New();
//	myCallback->Plane = plane;
//
//	//添加actor
//
//	renderer->AddActor(actor);
//	renderer->AddActor(cut_actor);
//
//	vtkSmartPointer<vtkImplicitPlaneRepresentation> rep =
//		vtkSmartPointer<vtkImplicitPlaneRepresentation>::New();
//	rep->SetPlaceFactor(1.25); // This must be set prior to placing the widget
//	rep->PlaceWidget(actor->GetBounds());
//	rep->SetNormal(plane->GetNormal());
//	rep->SetOrigin(plane->GetOrigin());
//	rep->OutlineTranslationOff();//锁定最外层边框
//
//	
//	vtkSmartPointer<vtkImplicitPlaneWidget2> planeWidget =
//		vtkSmartPointer<vtkImplicitPlaneWidget2>::New();
//	planeWidget->SetInteractor(renderWindowInteractor);
//	planeWidget->SetRepresentation(rep);
//	planeWidget->AddObserver(vtkCommand::InteractionEvent, myCallback);
//
//
//	planeWidget->On();
//	renderWindowInteractor->Initialize();
//	renderWindow->Render();
//
//
//	/*
//	// Setup a visualization pipeline
//	vtkSmartPointer<vtkPlane> plane =
//		vtkSmartPointer<vtkPlane>::New();
//	plane->SetOrigin(reader->GetOutput()->GetCenter());
//
//
//	auto clipper = vtkSmartPointer<vtkClipDataSet>::New();
//
//	clipper->SetClipFunction(plane);
//	clipper->InsideOutOn();
//	clipper->SetInputData(reader->GetOutput());
//	//clipper->SetValue(0.0);
//	clipper->GenerateClippedOutputOn();
//	clipper->Update();
//
//	vtkNew<vtkLookupTable> lut1;
//	lut1->SetHueRange(0.4, 0.9);
//	// Create a mapper and actor
//	//vtkSmartPointer<vtkPolyDataMapper> mapper =
//	//	vtkSmartPointer<vtkPolyDataMapper>::New();
//	vtkSmartPointer<vtkDataSetMapper>mapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	mapper->SetInputConnection(clipper->GetOutputPort());
//	mapper->SetScalarRange(reader->GetOutput()->GetScalarRange());
//	mapper->SetLookupTable(lut1);
//	mapper->SetColorModeToMapScalars();
//	vtkSmartPointer<vtkActor> actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	// A renderer and render window
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	renderer->AddActor(actor);
//
//	// An interactor
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderWindow->Render();
//
//	// The callback will do the work
//	vtkSmartPointer<vtkIPWCallback> myCallback =
//		vtkSmartPointer<vtkIPWCallback>::New();
//	myCallback->Plane = plane;
//	myCallback->Actor = actor;
//
//	vtkSmartPointer<vtkImplicitPlaneRepresentation> rep =
//		vtkSmartPointer<vtkImplicitPlaneRepresentation>::New();
//	rep->SetPlaceFactor(1.25); // This must be set prior to placing the widget
//	rep->PlaceWidget(actor->GetBounds());
//	rep->SetNormal(plane->GetNormal());
//	rep->SetOrigin(plane->GetOrigin());
//	rep->OutlineTranslationOff();//锁定最外层边框
//
//	vtkSmartPointer<vtkImplicitPlaneWidget2> planeWidget =
//		vtkSmartPointer<vtkImplicitPlaneWidget2>::New();
//	planeWidget->SetInteractor(renderWindowInteractor);
//	planeWidget->SetRepresentation(rep);
//	planeWidget->AddObserver(vtkCommand::InteractionEvent, myCallback);
//
//	// Render
//
//	renderWindowInteractor->Initialize();
//	renderWindow->Render();
//	planeWidget->On();
//
//	// Begin mouse interaction
//	renderWindowInteractor->Start();
//	*/
//
//	return EXIT_SUCCESS;
//}


//#include <vtkSmartPointer.h>
//
//#include <vtkXMLPolyDataReader.h>
//#include <vtkSphereSource.h>
//#include <vtkClipPolyData.h>
//#include <vtkPlane.h>
//
//#include <vtkCommand.h>
//#include <vtkImplicitPlaneWidget2.h>
//#include <vtkImplicitPlaneRepresentation.h>
//#include<vtkPolyDataMapper.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include<vtkUnstructuredGrid.h>
//#include<vtkUnstructuredGridReader.h>
//#include<vtkDataSetSurfaceFilter.h>
//#include<vtkDataSetMapper.h>
//#include<vtkNew.h>
//#include<vtkLookupTable.h>
//
//#include<vtkClipDataSet.h>
//#include<vtkCutter.h>
//
//int main(int argc, char* argv[])
//{
//
//	vtkSmartPointer<vtkUnstructuredGridReader> reader =
//		vtkSmartPointer<vtkUnstructuredGridReader>::New();
//	reader->SetFileName("third02.vtk");
//	reader->Update();
//
//	// Setup a visualization pipeline
//	vtkSmartPointer<vtkPlane> plane =
//		vtkSmartPointer<vtkPlane>::New();
//	plane->SetOrigin(reader->GetOutput()->GetCenter());
//	plane->SetNormal(1, 0, 0);
//
//	auto cutter = vtkSmartPointer<vtkCutter>::New();
//
//	cutter->SetInputConnection(reader->GetOutputPort());
//	cutter->SetCutFunction(plane);
//	
//
//	vtkNew<vtkLookupTable> lut1;
//	lut1->SetHueRange(0.4, 0.9);
//	// Create a mapper and actor
//	//vtkSmartPointer<vtkPolyDataMapper> mapper =
//	//	vtkSmartPointer<vtkPolyDataMapper>::New();
//	vtkSmartPointer<vtkDataSetMapper>mapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	mapper->SetInputConnection(cutter->GetOutputPort());
//	mapper->SetScalarRange(reader->GetOutput()->GetScalarRange());
//	mapper->SetLookupTable(lut1);
//	mapper->SetColorModeToMapScalars();
//	vtkSmartPointer<vtkActor> actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	// A renderer and render window
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	renderer->AddActor(actor);
//
//	// An interactor
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderWindow->Render();
//
//
//	// Render
//
//	renderWindowInteractor->Initialize();
//	renderWindow->Render();
//
//	// Begin mouse interaction
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}

//面单元的立方体

//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkCellArray.h>
//#include <vtkFloatArray.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPointData.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//
//#include <array>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	std::array<std::array<double, 3>, 8> pts = { {{{0, 0, 0}},
//												 {{1, 0, 0}},
//												 {{1, 1, 0}},
//												 {{0, 1, 0}},
//												 {{0, 0, 1}},
//												 {{1, 0, 1}},
//												 {{1, 1, 1}},
//												 {{0, 1, 1}}} };
//	// The ordering of the corner points on each face.
//	std::array<std::array<vtkIdType, 4>, 6> ordering = { {{{0, 1, 2, 3}},
//														 {{4, 5, 6, 7}},
//														 {{0, 1, 5, 4}},
//														 {{1, 2, 6, 5}},
//														 {{2, 3, 7, 6}},
//														 {{3, 0, 4, 7}}} };
//
//	// We'll create the building blocks of polydata including data attributes.
//	vtkNew<vtkPolyData> cube;
//	vtkNew<vtkPoints> points;
//	vtkNew<vtkCellArray> polys;
//	vtkNew<vtkFloatArray> scalars;
//
//	// Load the point, cell, and data attributes.
//	for (auto i = 0ul; i < pts.size(); ++i)
//	{
//		points->InsertPoint(i, pts[i].data());
//		scalars->InsertTuple1(i, i);
//	}
//	for (auto&& i : ordering)
//	{
//		polys->InsertNextCell(vtkIdType(i.size()), i.data());
//	}
//
//	// We now assign the pieces to the vtkPolyData.
//	cube->SetPoints(points);
//	cube->SetPolys(polys);
//	cube->GetPointData()->SetScalars(scalars);
//
//	// Now we'll look at it.
//	vtkNew<vtkPolyDataMapper> cubeMapper;
//	cubeMapper->SetInputData(cube);
//	cubeMapper->SetScalarRange(cube->GetScalarRange());
//	vtkNew<vtkActor> cubeActor;
//	cubeActor->SetMapper(cubeMapper);
//
//	// The usual rendering stuff.
//	vtkNew<vtkCamera> camera;
//	camera->SetPosition(1, 1, 1);
//	camera->SetFocalPoint(0, 0, 0);
//
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renWin;
//	renWin->AddRenderer(renderer);
//	renWin->SetWindowName("Cube");
//
//	vtkNew<vtkRenderWindowInteractor> iren;
//	iren->SetRenderWindow(renWin);
//
//	renderer->AddActor(cubeActor);
//	renderer->SetActiveCamera(camera);
//	renderer->ResetCamera();
//	renderer->SetBackground(colors->GetColor3d("Cornsilk").GetData());
//
//	renWin->SetSize(600, 600);
//
//	// interact with data
//	renWin->Render();
//	iren->Start();
//
//	return EXIT_SUCCESS;
//}


//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkPropCollection.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkInteractorStyleTrackball.h>

//int main(int, char* [])
//{
//	// Sphere 1
//	vtkSmartPointer<vtkSphereSource> sphereSource1 =
//		vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource1->SetCenter(0.0, 0.0, 0.0);
//	sphereSource1->SetRadius(4.0);
//	sphereSource1->Update();
//
//	vtkSmartPointer<vtkPolyDataMapper> mapper1 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper1->SetInputConnection(sphereSource1->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> actor1 =
//		vtkSmartPointer<vtkActor>::New();
//	actor1->SetMapper(mapper1);
//
//	// Sphere 2
//	vtkSmartPointer<vtkSphereSource> sphereSource2 =
//		vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource2->SetCenter(10.0, 0.0, 0.0);
//	sphereSource2->SetRadius(3.0);
//	sphereSource2->Update();
//
//	// Create a mapper
//	vtkSmartPointer<vtkPolyDataMapper> mapper2 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper2->SetInputConnection(sphereSource2->GetOutputPort());
//
//	// Create an actor
//	vtkSmartPointer<vtkActor> actor2 =
//		vtkSmartPointer<vtkActor>::New();
//	actor2->SetMapper(mapper2);
//
//	// A renderer and render window
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	// An interactor
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the scene
//	renderer->AddActor(actor1);
//	renderer->AddActor(actor2);
//	renderer->SetBackground(1, 1, 1); // Background color white
//
//	// Render an image (lights and cameras are created automatically)
//	renderWindow->Render();
//
//	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
//		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderWindowInteractor->Start();
//
//	// Set the background to red so we know we are done with the original two sphere display
//	renderer->SetBackground(1, 0, 0);
//
//	// Hide one actor at a time
//	vtkPropCollection* props = renderer->GetViewProps(); //iterate through and set each visibility to 0
//	props->InitTraversal();
//	for (int i = 0; i < props->GetNumberOfItems(); i++)
//	{
//		props->GetNextProp()->VisibilityOff();
//
//		renderer->ResetCamera();
//		renderWindow->Render();
//
//		renderWindowInteractor->Start();
//	}
//
//	return EXIT_SUCCESS;
//}


//手动输入点、单元显示，带色标、标量数据，贴吧https://tieba.baidu.com/p/2264333210?red_tag=2049147286

//#include "vtkCamera.h"
//#include "vtkCellDataToPointData.h"
//#include "vtkCompositeDataPipeline.h"
//#include "vtkContourFilter.h"
//#include "vtkDebugLeaks.h"
//#include "vtkHierarchicalDataSetGeometryFilter.h"
//#include "vtkOutlineCornerFilter.h"
//#include "vtkHierarchicalPolyDataMapper.h"
//#include "vtkProperty.h"
//#include "vtkRenderer.h"
//#include "vtkRenderWindow.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkFloatArray.h"
//#include "vtkCellArray.h"
//#include "vtkPointData.h"
//#include "vtkPolyDataMapper.h"
//#include "vtkIdFilter.h"
//#include "vtkBandedPolyDataContourFilter.h"
//#include "vtkLabeledDataMapper.h"
//#include "vtkActor2D.h" 
//int main(int argc, char* argv[])
//{
//	//# Manually create cells of various types: vertex, polyvertex, line,
//	//# polyline, triangle, quad, pentagon, and triangle strip.
//	vtkPoints* pts = vtkPoints::New();
//	pts->InsertPoint(0, 0, 0, 0);
//	pts->InsertPoint(1, 0, 1, 0);
//	pts->InsertPoint(2, 0, 2, 0);
//	pts->InsertPoint(3, 1, 0, 0);
//	pts->InsertPoint(4, 1, 1, 0);
//	pts->InsertPoint(5, 1, 2, 0);
//	pts->InsertPoint(6, 2, 0, 0);
//	pts->InsertPoint(7, 2, 2, 0);
//	pts->InsertPoint(8, 3, 0, 0);
//	pts->InsertPoint(9, 3, 1, 0);
//	pts->InsertPoint(10, 3, 2, 0);
//	pts->InsertPoint(11, 4, 0, 0);
//	pts->InsertPoint(12, 6, 0, 0);
//	pts->InsertPoint(13, 5, 2, 0);
//	pts->InsertPoint(14, 7, 0, 0);
//	pts->InsertPoint(15, 9, 0, 0);
//	pts->InsertPoint(16, 7, 2, 0);
//	pts->InsertPoint(17, 9, 2, 0);
//	pts->InsertPoint(18, 10, 0, 0);
//	pts->InsertPoint(19, 12, 0, 0);
//	pts->InsertPoint(20, 10, 1, 0);
//	pts->InsertPoint(21, 12, 1, 0);
//	pts->InsertPoint(22, 10, 2, 0);
//	pts->InsertPoint(23, 12, 2, 0);
//	pts->InsertPoint(24, 10, 3, 0);
//	pts->InsertPoint(25, 12, 3, 0);
//	pts->InsertPoint(26, 14, 0, 0);
//	pts->InsertPoint(27, 14, 2, 0);
//	pts->InsertPoint(28, 16, 2, 0);
//	pts->InsertPoint(29, 16, 0, 0);
//	pts->InsertPoint(30, 14, 0, 2);
//	pts->InsertPoint(31, 14, 2, 2);
//	pts->InsertPoint(32, 16, 2, 2);
//	pts->InsertPoint(33, 16, 0, 2); vtkCellArray* verts = vtkCellArray::New();
//	verts->InsertNextCell(1);
//	verts->InsertCellPoint(0);
//	verts->InsertNextCell(1);
//	verts->InsertCellPoint(1);
//	verts->InsertNextCell(1);
//	verts->InsertCellPoint(2);
//	verts->InsertNextCell(3);
//	verts->InsertCellPoint(3);
//	verts->InsertCellPoint(4);
//	verts->InsertCellPoint(5); vtkCellArray* lines = vtkCellArray::New();
//	lines->InsertNextCell(2);
//	lines->InsertCellPoint(6);
//	lines->InsertCellPoint(7);
//	lines->InsertNextCell(3);
//	lines->InsertCellPoint(8);
//	lines->InsertCellPoint(9);
//	lines->InsertCellPoint(10); vtkCellArray* polys = vtkCellArray::New();
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(14);
//	polys->InsertCellPoint(15);
//	polys->InsertCellPoint(17);
//	polys->InsertCellPoint(16);
//	polys->InsertNextCell(3);
//	polys->InsertCellPoint(11);
//	polys->InsertCellPoint(12);
//	polys->InsertCellPoint(13);
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(26);
//	polys->InsertCellPoint(27);
//	polys->InsertCellPoint(28);
//	polys->InsertCellPoint(29);
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(30);
//	polys->InsertCellPoint(31);
//	polys->InsertCellPoint(32);
//	polys->InsertCellPoint(33);
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(26);
//	polys->InsertCellPoint(27);
//	polys->InsertCellPoint(31);
//	polys->InsertCellPoint(30);
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(28);
//	polys->InsertCellPoint(29);
//	polys->InsertCellPoint(33);
//	polys->InsertCellPoint(32);
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(27);
//	polys->InsertCellPoint(28);
//	polys->InsertCellPoint(32);
//	polys->InsertCellPoint(31);
//	polys->InsertNextCell(4);
//	polys->InsertCellPoint(26);
//	polys->InsertCellPoint(29);
//	polys->InsertCellPoint(33);
//	polys->InsertCellPoint(30); vtkCellArray* strips = vtkCellArray::New();
//	strips->InsertNextCell(8);
//	strips->InsertCellPoint(19);
//	strips->InsertCellPoint(18);
//	strips->InsertCellPoint(21);
//	strips->InsertCellPoint(20);
//	strips->InsertCellPoint(23);
//	strips->InsertCellPoint(22);
//	strips->InsertCellPoint(25);
//	strips->InsertCellPoint(24); vtkFloatArray* scalars = vtkFloatArray::New();
//	scalars->SetNumberOfTuples(33);
//	scalars->SetTuple1(0, 0);
//	scalars->SetTuple1(1, 50);
//	scalars->SetTuple1(2, 100);
//	scalars->SetTuple1(3, 0);
//	scalars->SetTuple1(4, 50);
//	scalars->SetTuple1(5, 100);
//	scalars->SetTuple1(6, 10);
//	scalars->SetTuple1(7, 90);
//	scalars->SetTuple1(8, 10);
//	scalars->SetTuple1(9, 50);
//	scalars->SetTuple1(10, 90);
//	scalars->SetTuple1(11, 10);
//	scalars->SetTuple1(12, 40);
//	scalars->SetTuple1(13, 100);
//	scalars->SetTuple1(14, 0);
//	scalars->SetTuple1(15, 60);
//	scalars->SetTuple1(16, 40);
//	scalars->SetTuple1(17, 100);
//	scalars->SetTuple1(18, 0);
//	scalars->SetTuple1(19, 25);
//	scalars->SetTuple1(20, 25);
//	scalars->SetTuple1(21, 50);
//	scalars->SetTuple1(22, 50);
//	scalars->SetTuple1(23, 75);
//	scalars->SetTuple1(24, 75);
//	scalars->SetTuple1(25, 100);
//	scalars->SetTuple1(26, 0);
//	scalars->SetTuple1(27, 05);
//	scalars->SetTuple1(28, 25);
//	scalars->SetTuple1(29, 50);
//	scalars->SetTuple1(30, 50);
//	scalars->SetTuple1(31, 75);
//	scalars->SetTuple1(32, 75);
//	scalars->SetTuple1(33, 100); vtkPolyData* polyData = vtkPolyData::New();
//	polyData->SetPoints(pts);
//	polyData->SetVerts(verts);
//	polyData->SetLines(lines);
//	polyData->SetPolys(polys);
//	polyData->SetStrips(strips);
//	(polyData->GetPointData())->SetScalars(scalars); vtkBandedPolyDataContourFilter* bf = vtkBandedPolyDataContourFilter::New();
//	bf->SetInputData(polyData);
//	//bf->GenerateValues(3,25,75 );
//	//bf->GenerateValues(9,5,95 );
//	bf->SetNumberOfContours(9);
//	bf->SetValue(0, 5);
//	bf->SetValue(1, 15);
//	bf->SetValue(2, 25);
//	bf->SetValue(3, 35);
//	bf->SetValue(4, 45);
//	bf->SetValue(5, 55);
//	bf->SetValue(6, 65);
//	bf->SetValue(7, 75);
//	bf->SetValue(8, 85);
//	vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
//	mapper->SetInputConnection(bf->GetOutputPort());
//	mapper->SetScalarModeToUseCellData();
//	//mapper->SetScalarRange(0,4);
//	mapper->SetScalarRange(0, 10);
//	vtkActor* actor = vtkActor::New();
//	actor->SetMapper(mapper); vtkIdFilter* ids = vtkIdFilter::New();
//	ids->SetInputConnection(bf->GetOutputPort());
//	ids->PointIdsOn();
//	ids->CellIdsOn();
//	ids->FieldDataOn();
//	vtkLabeledDataMapper* ldm = vtkLabeledDataMapper::New();
//	ldm->SetInputConnection(ids->GetOutputPort());
//	//# ldm SetLabelFormat "%g"
//	ldm->SetLabelModeToLabelFieldData();
//	vtkActor2D* pointLabels = vtkActor2D::New();
//	pointLabels->SetMapper(ldm);
//	//# Create the RenderWindow, Renderer and both Actors
//	//#
//	vtkRenderer* ren1 = vtkRenderer::New();
//	vtkRenderWindow* renWin = vtkRenderWindow::New();
//	renWin->AddRenderer(ren1);
//	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
//	iren->SetRenderWindow(renWin); //# Add the actors to the renderer, set the background and size
//	//#
//	ren1->AddActor(actor);
//	ren1->AddActor2D(pointLabels);
//	//#for debugging only
//	ren1->SetBackground(0, 0, 0);
//	renWin->SetSize(300, 80);
//	renWin->Render();
//	(ren1->GetActiveCamera())->Zoom(3);
//	renWin->Render();
//	iren->Initialize();
//	iren->Start();
//	return 0;
//}


//选中actor，加外边框

//#include <vtkAutoInit.h>
//#include <vtkObject.h>
//#include <vtkSmartPointer.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkPolyDataSilhouette.h>
//#include <vtkMath.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPropPicker.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSphereSource.h>
//#include <vtkNamedColors.h>
//
//namespace
//{
//	// Handle mouse events
//	class MouseInteractorHighLightActor : public vtkInteractorStyleTrackballCamera
//	{
//	public:
//		static MouseInteractorHighLightActor* New();
//		vtkTypeMacro(MouseInteractorHighLightActor,
//			vtkInteractorStyleTrackballCamera);
//
//		MouseInteractorHighLightActor()
//		{
//			LastPickedActor = nullptr;
//			SilhouetteActor = nullptr;
//			Silhouette = nullptr;
//		}
//		virtual ~MouseInteractorHighLightActor()
//		{
//		}
//		virtual void OnLeftButtonDown() override
//		{
//			int* clickPos = this->GetInteractor()->GetEventPosition();
//
//			// Pick from this location.
//			auto picker =
//				vtkSmartPointer<vtkPropPicker>::New();
//			picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
//			//this->LastPickedActor = picker->GetActor();
//
//			// If we picked something before, remove the silhouette actor and
//			// generate a new one
//			if (this->LastPickedActor)
//			{
//				this->GetDefaultRenderer()->RemoveActor(this->SilhouetteActor);
//			}
//			this->LastPickedActor = picker->GetActor();
//			if (this->LastPickedActor)
//			{
//				this->GetDefaultRenderer()->RemoveActor(this->SilhouetteActor);
//				
//				// Highlight the picked actor by generating a silouhette
//				this->Silhouette->SetInputData(
//					dynamic_cast<vtkPolyDataMapper*>(this->LastPickedActor->GetMapper())->GetInput());
//				this->GetDefaultRenderer()->AddActor(this->SilhouetteActor);
//			}
//
//			// Forward events
//			vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
//		}
//		void SetSilhouette(vtkPolyDataSilhouette* silhouette)
//		{
//			this->Silhouette = silhouette;
//		}
//		void SetSilhouetteActor(vtkActor* silhouetteActor)
//		{
//			this->SilhouetteActor = silhouetteActor;
//		}
//	private:
//		vtkActor* LastPickedActor;
//		vtkActor* SilhouetteActor;
//		vtkPolyDataSilhouette* Silhouette;
//	};
//}
//
//vtkStandardNewMacro(MouseInteractorHighLightActor);
//
//// Execute application.
//int main(int argc, char* argv[])
//{
//	auto colors =
//		vtkSmartPointer<vtkNamedColors>::New();
//	colors->SetColor("Bkg", 0.3, 0.4, 0.5);
//
//	int numberOfSpheres = 10;
//	if (argc > 1)
//	{
//		numberOfSpheres = atoi(argv[1]);
//	}
//	// Create a renderer and render window
//	auto renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	auto renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	// An interactor
//	auto interactor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	interactor->SetRenderWindow(renderWindow);
//
//	for (int i = 0; i < numberOfSpheres; ++i)
//	{
//		auto source =
//			vtkSmartPointer<vtkSphereSource>::New();
//		double x, y, z, radius;
//		x = vtkMath::Random(-5, 5);
//		y = vtkMath::Random(-5, 5);
//		z = vtkMath::Random(-5, 5);
//		radius = vtkMath::Random(0.5, 1.0);
//		source->SetRadius(radius);
//		source->SetCenter(x, y, z);
//		source->SetPhiResolution(11);
//		source->SetThetaResolution(21);
//		auto mapper =
//			vtkSmartPointer<vtkPolyDataMapper>::New();
//		mapper->SetInputConnection(source->GetOutputPort());
//		auto actor =
//			vtkSmartPointer<vtkActor>::New();
//		actor->SetMapper(mapper);
//		double r, g, b;
//		r = vtkMath::Random(0.4, 1.0);
//		g = vtkMath::Random(0.4, 1.0);
//		b = vtkMath::Random(0.4, 1.0);
//		actor->GetProperty()->SetDiffuseColor(r, g, b);
//		actor->GetProperty()->SetDiffuse(0.8);
//		actor->GetProperty()->SetSpecular(0.5);
//		actor->GetProperty()->SetSpecularColor(
//			colors->GetColor3d("White").GetData());
//		actor->GetProperty()->SetSpecularPower(30.0);
//		renderer->AddActor(actor);
//	}
//
//	renderer->SetBackground(colors->GetColor3d("Bkg").GetData());
//
//	// Render and interact
//	renderWindow->Render();
//
//	// Create the silhouette pipeline, the input data will be set in the
//	// interactor
//	auto silhouette =
//		vtkSmartPointer<vtkPolyDataSilhouette>::New();
//	silhouette->SetCamera(renderer->GetActiveCamera());
//
//	// Create mapper and actor for silhouette
//	auto silhouetteMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	silhouetteMapper->SetInputConnection(silhouette->GetOutputPort());
//
//	auto silhouetteActor =
//		vtkSmartPointer<vtkActor>::New();
//	silhouetteActor->SetMapper(silhouetteMapper);
//	silhouetteActor->GetProperty()->SetColor(colors->GetColor3d("Tomato").GetData());
//	silhouetteActor->GetProperty()->SetLineWidth(5);
//
//	// Set the custom type to use for interaction.
//	auto style =
//		vtkSmartPointer<MouseInteractorHighLightActor>::New();
//	style->SetDefaultRenderer(renderer);
//	style->SetSilhouetteActor(silhouetteActor);
//	style->SetSilhouette(silhouette);
//
//	interactor->SetInteractorStyle(style);
//
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}


//#include <vtkAutoInit.h>
//#include <vtkObject.h>
//#include <vtkSmartPointer.h>
//#include <vtkActor.h>
//#include <vtkCellArray.h>
//#include <vtkInteractorStyleTrackballActor.h>
//#include <vtkObjectFactory.h>
//#include <vtkCubeSource.h>
//#include <vtkSphereSource.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPropPicker.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//
//// Handle mouse events
//class MouseInteractorStyle5 : public vtkInteractorStyleTrackballActor
//{
//public:
//	static MouseInteractorStyle5* New();
//	vtkTypeMacro(MouseInteractorStyle5, vtkInteractorStyleTrackballActor);
//
//	virtual void OnLeftButtonDown() override
//	{
//		// Forward events
//		vtkInteractorStyleTrackballActor::OnLeftButtonDown();
//
//		this->Cube->GetProperty()->SetColor(.5, .5, 1);
//		this->Sphere->GetProperty()->SetColor(1, 1, 1);
//		if (this->InteractionProp == this->Cube)
//		{
//			std::cout << "Picked cube." << std::endl;
//			this->Cube->GetProperty()->SetColor(1, 0, 0);
//		}
//		else if (this->InteractionProp == this->Sphere)
//		{
//			std::cout << "Picked sphere." << std::endl;
//			this->Sphere->GetProperty()->SetColor(0, 1, 0);
//		}
//	}
//
//	vtkActor* Cube;
//	vtkActor* Sphere;
//};
//
//vtkStandardNewMacro(MouseInteractorStyle5);
//
//int main(int, char* [])
//{
//	// Create a cube
//	vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
//	cubeSource->Update();
//
//	vtkSmartPointer<vtkPolyDataMapper> cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> cubeActor = vtkSmartPointer<vtkActor>::New();
//	cubeActor->SetMapper(cubeMapper);
//
//	// Create a sphere
//	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource->SetCenter(5, 0, 0);
//	sphereSource->Update();
//
//	// Create a mapper
//	vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
//
//	// Create an actor
//	vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
//	sphereActor->SetMapper(sphereMapper);
//
//	// A renderer and render window
//	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	// An interactor
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Set the custom stype to use for interaction.
//	vtkSmartPointer<MouseInteractorStyle5> style = vtkSmartPointer<MouseInteractorStyle5>::New();
//	style->SetDefaultRenderer(renderer);
//	style->Cube = cubeActor;
//	style->Sphere = sphereActor;
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderer->AddActor(cubeActor);
//	renderer->AddActor(sphereActor);
//	renderer->SetBackground(.1, .2, .3);
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Initialize();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//只控制单个actor，不控制renderer

//#include <vtkActor.h>
//#include <vtkInteractorStyleTrackballActor.h>
//#include <vtkNew.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkNew.h>
//#include <vtkSphereSource.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	// Sphere 1
//	vtkNew<vtkSphereSource> sphereSource1;
//	sphereSource1->SetCenter(0.0, 0.0, 0.0);
//	sphereSource1->SetRadius(4.0);
//	sphereSource1->Update();
//
//	vtkNew<vtkPolyDataMapper> mapper1;
//	mapper1->SetInputConnection(sphereSource1->GetOutputPort());
//
//	vtkNew<vtkActor> actor1;
//	actor1->SetMapper(mapper1);
//	actor1->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
//
//	// Sphere 2
//	vtkNew<vtkSphereSource> sphereSource2;
//	sphereSource2->SetCenter(10.0, 0.0, 0.0);
//	sphereSource2->SetRadius(3.0);
//	sphereSource2->Update();
//
//	// Create a mapper
//	vtkNew<vtkPolyDataMapper> mapper2;
//	mapper2->SetInputConnection(sphereSource2->GetOutputPort());
//
//	// Create an actor
//	vtkNew<vtkActor> actor2;
//	actor2->SetMapper(mapper2);
//	actor2->GetProperty()->SetColor(colors->GetColor3d("Cornsilk").GetData());
//
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("MoveActor");
//
//	// An interactor
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the scene
//	renderer->AddActor(actor1);
//	renderer->AddActor(actor2);
//	renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());
//
//	// Render an image (lights and cameras are created automatically)
//	renderWindow->Render();
//
//	vtkNew<vtkInteractorStyleTrackballActor> style;
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	// Begin mouse interaction
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//矩形拾取点

//#include <vtkActor.h>
//#include <vtkAreaPicker.h>
//#include <vtkDataSetMapper.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkExtractGeometry.h>
//#include <vtkIdFilter.h>
//#include <vtkIdTypeArray.h>
//#include <vtkInteractorStyleRubberBandPick.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkObjectFactory.h>
//#include <vtkPlanes.h>
//#include <vtkPointData.h>
//#include <vtkPointSource.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkRendererCollection.h>
//#include <vtkSmartPointer.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkVersion.h>
//#include <vtkVertexGlyphFilter.h>
//
//
//#if VTK_VERSION_NUMBER >= 89000000000ULL
//#define VTK890 1
//#endif
//
//namespace {
//	// Define interaction style
//	class InteractorStyle : public vtkInteractorStyleRubberBandPick
//	{
//	public:
//		static InteractorStyle* New();
//		vtkTypeMacro(InteractorStyle, vtkInteractorStyleRubberBandPick);
//
//		InteractorStyle()
//		{
//			this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
//			this->SelectedActor = vtkSmartPointer<vtkActor>::New();
//			this->SelectedActor->SetMapper(SelectedMapper);
//		}
//
//		virtual void OnLeftButtonUp() override
//		{
//			vtkNew<vtkNamedColors> colors;
//
//			// Forward events
//			vtkInteractorStyleRubberBandPick::OnLeftButtonUp();
//
//			vtkPlanes* frustum =
//				static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())
//				->GetFrustum();
//
//			vtkNew<vtkExtractGeometry> extractGeometry;
//			extractGeometry->SetImplicitFunction(frustum);
//			extractGeometry->SetInputData(this->Points);
//			extractGeometry->Update();
//
//			vtkNew<vtkVertexGlyphFilter> glyphFilter;
//			glyphFilter->SetInputConnection(extractGeometry->GetOutputPort());
//			glyphFilter->Update();
//
//			vtkPolyData* selected = glyphFilter->GetOutput();
//			std::cout << "Selected " << selected->GetNumberOfPoints() << " points."
//				<< std::endl;
//			std::cout << "Selected " << selected->GetNumberOfCells() << " cells."
//				<< std::endl;
//			this->SelectedMapper->SetInputData(selected);
//			this->SelectedMapper->ScalarVisibilityOff();
//
//			vtkIdTypeArray* ids = dynamic_cast<vtkIdTypeArray*>(
//				selected->GetPointData()->GetArray("OriginalIds"));
//			for (vtkIdType i = 0; i < ids->GetNumberOfTuples(); i++)
//			{
//				std::cout << "Id " << i << " : " << ids->GetValue(i) << std::endl;
//			}
//
//			this->SelectedActor->GetProperty()->SetColor(
//				colors->GetColor3d("Red").GetData());
//			this->SelectedActor->GetProperty()->SetPointSize(5);
//
//			this->CurrentRenderer->AddActor(SelectedActor);
//			this->GetInteractor()->GetRenderWindow()->Render();
//			this->HighlightProp(NULL);
//		}
//
//		void SetPoints(vtkSmartPointer<vtkPolyData> points)
//		{
//			this->Points = points;
//		}
//
//	private:
//		vtkSmartPointer<vtkPolyData> Points;
//		vtkSmartPointer<vtkActor> SelectedActor;
//		vtkSmartPointer<vtkDataSetMapper> SelectedMapper;
//	};
//
//	vtkStandardNewMacro(InteractorStyle);
//} // namespace
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkPointSource> pointSource;
//	pointSource->SetNumberOfPoints(20);
//	pointSource->Update();
//
//	vtkNew<vtkIdFilter> idFilter;
//	idFilter->SetInputConnection(pointSource->GetOutputPort());
//#if VTK890
//	idFilter->SetCellIdsArrayName("OriginalIds");
//	idFilter->SetPointIdsArrayName("OriginalIds");
//#else
//	idFilter->SetIdsArrayName("OriginalIds");
//#endif
//	idFilter->Update();
//
//	vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
//	surfaceFilter->SetInputConnection(idFilter->GetOutputPort());
//	surfaceFilter->Update();
//
//	vtkPolyData* input = surfaceFilter->GetOutput();
//
//	// Create a mapper and actor
//	vtkNew<vtkPolyDataMapper> mapper;
//	mapper->SetInputData(input);
//	mapper->ScalarVisibilityOff();
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetPointSize(3);
//	actor->GetProperty()->SetColor(colors->GetColor3d("Gold").GetData());
//
//	// Visualize
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("HighlightSelectedPoints");
//
//	vtkNew<vtkAreaPicker> areaPicker;
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetPicker(areaPicker);
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(actor);
//	renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());
//
//	renderWindow->Render();
//
//	vtkNew<InteractorStyle> style;
//	style->SetPoints(input);
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//矩形拾取单元 类似procast框选，可以选择被遮挡部分，前处理
//Press 'r' to enter selection mode.
//HighlightSelection

//#include <vtkActor.h>
//#include <vtkAreaPicker.h>
//#include <vtkDataSetMapper.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkExtractPolyDataGeometry.h>
//#include <vtkIdFilter.h>
//#include <vtkIdTypeArray.h>
//#include <vtkInteractorStyleRubberBandPick.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkObjectFactory.h>
//#include <vtkPlanes.h>
//#include <vtkPointData.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkRendererCollection.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkVersion.h>
//#include <vtkVertexGlyphFilter.h>
//
//#if VTK_VERSION_NUMBER >= 89000000000ULL
//#define VTK890 1
//#endif
//
//#include <vtkBYUReader.h>
//#include <vtkOBJReader.h>
//#include <vtkPLYReader.h>
//#include <vtkPolyDataReader.h>
//#include <vtkSTLReader.h>
//#include <vtkXMLPolyDataReader.h>
//#include <vtksys/SystemTools.hxx>
//
//#define VTKISRBP_ORIENT 0
//#define VTKISRBP_SELECT 1
//
//namespace {
//	// Define interaction style
//	class HighlightInteractorStyle : public vtkInteractorStyleRubberBandPick
//	{
//	public:
//		static HighlightInteractorStyle* New();
//		vtkTypeMacro(HighlightInteractorStyle, vtkInteractorStyleRubberBandPick);
//
//		HighlightInteractorStyle() : vtkInteractorStyleRubberBandPick()
//		{
//			this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
//			this->SelectedActor = vtkSmartPointer<vtkActor>::New();
//			this->SelectedActor->SetMapper(SelectedMapper);
//		}
//
//		virtual void OnLeftButtonUp() override
//		{
//			// Forward events
//			vtkInteractorStyleRubberBandPick::OnLeftButtonUp();
//
//			if (this->CurrentMode == VTKISRBP_SELECT)
//			{
//				vtkNew<vtkNamedColors> colors;
//
//				vtkPlanes* frustum =
//					static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())
//					->GetFrustum();
//
//				vtkNew<vtkExtractPolyDataGeometry> extractPolyDataGeometry;
//				extractPolyDataGeometry->SetInputData(this->PolyData);
//				extractPolyDataGeometry->SetImplicitFunction(frustum);
//				extractPolyDataGeometry->Update();
//
//				std::cout << "Extracted "
//					<< extractPolyDataGeometry->GetOutput()->GetNumberOfCells()
//					<< " cells." << std::endl;
//				this->SelectedMapper->SetInputData(extractPolyDataGeometry->GetOutput());
//				this->SelectedMapper->ScalarVisibilityOff();
//
//				//        vtkIdTypeArray* ids =
//				//        dynamic_cast<vtkIdTypeArray*>(selected->GetPointData()->GetArray("OriginalIds"));
//
//				this->SelectedActor->GetProperty()->SetColor(
//					colors->GetColor3d("Tomato").GetData());
//				this->SelectedActor->GetProperty()->SetPointSize(5);
//				this->SelectedActor->GetProperty()->SetRepresentationToWireframe();
//
//				this->GetInteractor()
//					->GetRenderWindow()
//					->GetRenderers()
//					->GetFirstRenderer()
//					->AddActor(SelectedActor);
//				this->GetInteractor()->GetRenderWindow()->Render();
//				this->HighlightProp(NULL);
//			}
//		}
//
//		void SetPolyData(vtkSmartPointer<vtkPolyData> polyData)
//		{
//			this->PolyData = polyData;
//		}
//
//	private:
//		vtkSmartPointer<vtkPolyData> PolyData;
//		vtkSmartPointer<vtkActor> SelectedActor;
//		vtkSmartPointer<vtkDataSetMapper> SelectedMapper;
//	};
//	vtkStandardNewMacro(HighlightInteractorStyle);
//
//	vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName);
//} // namespace
//
//int main(int argc, char* argv[])
//{
//	auto polyData = ReadPolyData(argc > 1 ? argv[1] : "");
//
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkIdFilter> idFilter;
//	idFilter->SetInputData(polyData);
//#if VTK890
//	idFilter->SetCellIdsArrayName("OriginalIds");
//	idFilter->SetPointIdsArrayName("OriginalIds");
//#else
//	idFilter->SetIdsArrayName("OriginalIds");
//#endif
//	idFilter->Update();
//
//	// This is needed to convert the ouput of vtkIdFilter (vtkDataSet) back to
//	// vtkPolyData
//	vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
//	surfaceFilter->SetInputConnection(idFilter->GetOutputPort());
//	surfaceFilter->Update();
//
//	vtkPolyData* input = surfaceFilter->GetOutput();
//
//	// Create a mapper and actor
//	vtkNew<vtkPolyDataMapper> mapper;
//	mapper->SetInputData(polyData);
//	mapper->ScalarVisibilityOff();
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetPointSize(5);
//	actor->GetProperty()->SetDiffuseColor(
//		colors->GetColor3d("Peacock").GetData());
//	actor->GetProperty()->EdgeVisibilityOn();
//	// Visualize
//	vtkNew<vtkRenderer> renderer;
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetSize(640, 480);
//	renderWindow->SetWindowName("HighlightSelection");
//
//	vtkNew<vtkAreaPicker> areaPicker;
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetPicker(areaPicker);
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(actor);
//	renderer->SetBackground(colors->GetColor3d("Tan").GetData());
//
//	renderWindow->Render();
//
//	vtkNew<HighlightInteractorStyle> style;
//	style->SetPolyData(input);
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}
//namespace {
//	vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName)
//	{
//		vtkSmartPointer<vtkPolyData> polyData;
//		std::string extension =
//			vtksys::SystemTools::GetFilenameLastExtension(std::string(fileName));
//		if (extension == ".ply")
//		{
//			vtkNew<vtkPLYReader> reader;
//			reader->SetFileName(fileName);
//			reader->Update();
//			polyData = reader->GetOutput();
//		}
//		else if (extension == ".vtp")
//		{
//			vtkNew<vtkXMLPolyDataReader> reader;
//			reader->SetFileName(fileName);
//			reader->Update();
//			polyData = reader->GetOutput();
//		}
//		else if (extension == ".obj")
//		{
//			vtkNew<vtkOBJReader> reader;
//			reader->SetFileName(fileName);
//			reader->Update();
//			polyData = reader->GetOutput();
//		}
//		else if (extension == ".stl")
//		{
//			vtkNew<vtkSTLReader> reader;
//			reader->SetFileName(fileName);
//			reader->Update();
//			polyData = reader->GetOutput();
//		}
//		else if (extension == ".vtk")
//		{
//			vtkNew<vtkPolyDataReader> reader;
//			reader->SetFileName(fileName);
//			reader->Update();
//			polyData = reader->GetOutput();
//		}
//		else if (extension == ".g")
//		{
//			vtkNew<vtkBYUReader> reader;
//			reader->SetGeometryFileName(fileName);
//			reader->Update();
//			polyData = reader->GetOutput();
//		}
//		else
//		{
//			vtkNew<vtkSphereSource> source;
//			source->SetPhiResolution(21);
//			source->SetThetaResolution(40);
//			source->Update();
//			polyData = source->GetOutput();
//		}
//		return polyData;
//	}
//} // namespace



//鼠标点击放置球体

//#include <vtkactor.h>
//#include <vtkcellarray.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtknamedcolors.h>
//#include <vtknew.h>
//#include <vtkobjectfactory.h>
//#include <vtkplanesource.h>
//#include <vtkpoints.h>
//#include <vtkpolydata.h>
//#include <vtkpolydatamapper.h>
//#include <vtkproppicker.h>
//#include <vtkproperty.h>
//#include <vtkrenderwindow.h>
//#include <vtkrenderwindowinteractor.h>
//#include <vtkrenderer.h>
//#include <vtkrenderercollection.h>
//#include <vtkspheresource.h>
//
//namespace {
//
//	// handle mouse events
//	class mouseinteractorstyle2 : public vtkInteractorStyleTrackballCamera
//	{
//	public:
//		static mouseinteractorstyle2* New();
//		vtkTypeMacro(mouseinteractorstyle2, vtkInteractorStyleTrackballCamera);
//		vtkNew<vtkNamedColors> colors;
//
//		virtual void OnLeftButtonDown() override
//		{
//			int* clickpos = this->GetInteractor()->GetEventPosition();
//
//			// pick from this location.
//			vtkNew<vtkPropPicker> picker;
//			picker->Pick(clickpos[0], clickpos[1], 0, this->GetDefaultRenderer());
//
//			double* pos = picker->GetPickPosition();
//			std::cout << "pick position (world coordinates) is: " << pos[0] << " "
//				<< pos[1] << " " << pos[2] << std::endl;
//			 
//			auto pickedactor = picker->GetActor();
//			if (pickedactor == nullptr)
//			{
//				std::cout << "no actor picked." << std::endl;
//			}
//			else
//			{
//				std::cout << "picked actor: " << picker->GetActor() << std::endl;
//				// create a sphere
//				// create a sphere
//				vtkNew<vtkSphereSource> spheresource;
//				spheresource->SetCenter(pos[0], pos[1], pos[2]);
//				spheresource->SetRadius(0.1);
//
//				// create a mapper and actor
//				vtkNew<vtkPolyDataMapper> mapper;
//				mapper->SetInputConnection(spheresource->GetOutputPort());
//
//				vtkNew<vtkActor> actor;
//				actor->SetMapper(mapper);
//				actor->GetProperty()->SetColor(colors->GetColor3d("mistyrose").GetData());
//
//				this->GetDefaultRenderer()->AddActor(actor);
//				// forward events
//				vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
//			}
//		}
//
//	private:
//	};
//
//	vtkStandardNewMacro(mouseinteractorstyle2);
//
//} // namespace
//
//// execute application.
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkPlaneSource> planesource;
//	planesource->Update();
//
//	// create a polydata object
//	vtkPolyData* polydata = planesource->GetOutput();
//
//	// create a mapper
//	vtkNew<vtkPolyDataMapper> mapper;
//	mapper->SetInputData(polydata);
//
//	// create an actor
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetColor(
//		colors->GetColor3d("lightgoldenrodyellow").GetData());
//
//	std::cout << "actor address: " << actor << std::endl;
//
//	// a renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderwindow;
//	renderwindow->AddRenderer(renderer);
//	renderwindow->SetWindowName("picking");
//
//	// an interactor
//	vtkNew<vtkRenderWindowInteractor> renderwindowinteractor;
//	renderwindowinteractor->SetRenderWindow(renderwindow);
//
//	// set the custom stype to use for interaction.
//	vtkNew<mouseinteractorstyle2> style;
//	style->SetDefaultRenderer(renderer);
//
//	renderwindowinteractor->SetInteractorStyle(style);
//
//	// add the actors to the scene
//	renderer->AddActor(actor);
//	renderer->SetBackground(colors->GetColor3d("dodgerblue").GetData());
//
//	// render and interact
//	renderwindow->Render();
//	renderwindowinteractor->Initialize();
//	renderwindowinteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//选中actor高亮显示HighlightPickedActor

//#include <vtkActor.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkMinimalStandardRandomSequence.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkObjectFactory.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPropPicker.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSphereSource.h>
//
//namespace {
//	// Handle mouse events
//	class MouseInteractorHighLightActor : public vtkInteractorStyleTrackballCamera
//	{
//	public:
//		static MouseInteractorHighLightActor* New();
//		vtkTypeMacro(MouseInteractorHighLightActor,
//			vtkInteractorStyleTrackballCamera);
//
//		MouseInteractorHighLightActor()
//		{
//			LastPickedActor = NULL;
//			LastPickedProperty = vtkProperty::New();
//		}
//		virtual ~MouseInteractorHighLightActor()
//		{
//			LastPickedProperty->Delete();
//		}
//		virtual void OnLeftButtonDown() override
//		{
//			vtkNew<vtkNamedColors> colors;
//
//			int* clickPos = this->GetInteractor()->GetEventPosition();
//
//			// Pick from this location.
//			vtkNew<vtkPropPicker> picker;
//			picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
//
//			// If we picked something before, reset its property
//			if (this->LastPickedActor)
//			{
//				this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
//			}
//			this->LastPickedActor = picker->GetActor();
//			if (this->LastPickedActor)
//			{
//				// Save the property of the picked actor so that we can
//				// restore it next time
//				this->LastPickedProperty->DeepCopy(this->LastPickedActor->GetProperty());
//				// Highlight the picked actor by changing its properties
//				this->LastPickedActor->GetProperty()->SetColor(
//					colors->GetColor3d("Red").GetData());
//				this->LastPickedActor->GetProperty()->SetDiffuse(1.0);
//				this->LastPickedActor->GetProperty()->SetSpecular(0.0);
//				this->LastPickedActor->GetProperty()->EdgeVisibilityOn();
//			}
//
//			// Forward events
//			vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
//		}
//
//	private:
//		vtkActor* LastPickedActor;
//		vtkProperty* LastPickedProperty;
//	};
//
//	vtkStandardNewMacro(MouseInteractorHighLightActor);
//} // namespace
//
//// Execute application.
//int main(int argc, char* argv[])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	int numberOfSpheres = 10;
//	if (argc > 1)
//	{
//		numberOfSpheres = atoi(argv[1]);
//	}
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("HighlightPickedActor");
//
//	// An interactor
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Set the custom type to use for interaction.
//	vtkNew<MouseInteractorHighLightActor> style;
//	style->SetDefaultRenderer(renderer);
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	vtkNew<vtkMinimalStandardRandomSequence> randomSequence;
//	randomSequence->SetSeed(8775070);
//	for (int i = 0; i < numberOfSpheres; ++i)
//	{
//		vtkNew<vtkSphereSource> source;
//		double x, y, z, radius;
//		// random position and radius
//		x = randomSequence->GetRangeValue(-5.0, 5.0);
//		randomSequence->Next();
//		y = randomSequence->GetRangeValue(-5.0, 5.0);
//		randomSequence->Next();
//		z = randomSequence->GetRangeValue(-5.0, 5.0);
//		randomSequence->Next();
//		radius = randomSequence->GetRangeValue(0.5, 1.0);
//		randomSequence->Next();
//		source->SetRadius(radius);
//		source->SetCenter(x, y, z);
//		source->SetPhiResolution(11);
//		source->SetThetaResolution(21);
//		vtkNew<vtkPolyDataMapper> mapper;
//		mapper->SetInputConnection(source->GetOutputPort());
//		vtkNew<vtkActor> actor;
//		actor->SetMapper(mapper);
//		double r, g, b;
//		r = randomSequence->GetRangeValue(0.4, 1.0);
//		randomSequence->Next();
//		g = randomSequence->GetRangeValue(0.4, 1.0);
//		randomSequence->Next();
//		b = randomSequence->GetRangeValue(0.4, 1.0);
//		randomSequence->Next();
//		actor->GetProperty()->SetDiffuseColor(r, g, b);
//		actor->GetProperty()->SetDiffuse(0.8);
//		actor->GetProperty()->SetSpecular(0.5);
//		actor->GetProperty()->SetSpecularColor(
//			colors->GetColor3d("White").GetData());
//		actor->GetProperty()->SetSpecularPower(30.0);
//		renderer->AddActor(actor);
//	}
//
//	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Initialize();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//#include <vtkActor.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkMinimalStandardRandomSequence.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkObjectFactory.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPropPicker.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSphereSource.h>
//
//namespace {
//	// Handle mouse events
//	class MouseInteractorHighLightActor : public vtkInteractorStyleTrackballCamera
//	{
//	public:
//		static MouseInteractorHighLightActor* New();
//		vtkTypeMacro(MouseInteractorHighLightActor,
//			vtkInteractorStyleTrackballCamera);
//
//		MouseInteractorHighLightActor()
//		{
//			LastPickedActor = NULL;
//			LastPickedProperty = vtkProperty::New();
//		}
//		virtual ~MouseInteractorHighLightActor()
//		{
//			LastPickedProperty->Delete();
//		}
//		virtual void OnLeftButtonDown() override
//		{
//			vtkNew<vtkNamedColors> colors;
//
//			int* clickPos = this->GetInteractor()->GetEventPosition();
//
//			// Pick from this location.
//			vtkNew<vtkPropPicker> picker;
//			picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
//			this->LastPickedActor = picker->GetActor();
//			// If we picked something before, reset its property
//			if (this->LastPickedActor)
//			{
//				this->LastPickedActor->GetProperty()->DeepCopy(this->LastPickedProperty);
//			}
//
//			if (this->LastPickedActor)
//			{
//				// Save the property of the picked actor so that we can
//				// restore it next time
//				this->LastPickedProperty->DeepCopy(this->LastPickedActor->GetProperty());
//				// Highlight the picked actor by changing its properties
//
//				//this->LastPickedActor->GetProperty()->SetColor(
//				//	colors->GetColor3d("Red").GetData());
//				//this->LastPickedActor->GetProperty()->SetDiffuse(1.0);
//				//this->LastPickedActor->GetProperty()->SetSpecular(0.0);
//				//this->LastPickedActor->GetProperty()->EdgeVisibilityOn();
//				this->GetDefaultRenderer()->RemoveActor(this->LastPickedActor);
//			}
//
//			// Forward events
//			vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
//		}
//
//	private:
//		vtkActor* LastPickedActor;
//		vtkProperty* LastPickedProperty;
//	};
//
//	vtkStandardNewMacro(MouseInteractorHighLightActor);
//} // namespace
//
//// Execute application.
//int main(int argc, char* argv[])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	int numberOfSpheres = 10;
//	if (argc > 1)
//	{
//		numberOfSpheres = atoi(argv[1]);
//	}
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("HighlightPickedActor");
//
//	// An interactor
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Set the custom type to use for interaction.
//	vtkNew<MouseInteractorHighLightActor> style;
//	style->SetDefaultRenderer(renderer);
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	vtkNew<vtkMinimalStandardRandomSequence> randomSequence;
//	randomSequence->SetSeed(8775070);
//	for (int i = 0; i < numberOfSpheres; ++i)
//	{
//		vtkNew<vtkSphereSource> source;
//		double x, y, z, radius;
//		// random position and radius
//		x = randomSequence->GetRangeValue(-5.0, 5.0);
//		randomSequence->Next();
//		y = randomSequence->GetRangeValue(-5.0, 5.0);
//		randomSequence->Next();
//		z = randomSequence->GetRangeValue(-5.0, 5.0);
//		randomSequence->Next();
//		radius = randomSequence->GetRangeValue(0.5, 1.0);
//		randomSequence->Next();
//		source->SetRadius(radius);
//		source->SetCenter(x, y, z);
//		source->SetPhiResolution(11);
//		source->SetThetaResolution(21);
//		vtkNew<vtkPolyDataMapper> mapper;
//		mapper->SetInputConnection(source->GetOutputPort());
//		vtkNew<vtkActor> actor;
//		actor->SetMapper(mapper);
//		double r, g, b;
//		r = randomSequence->GetRangeValue(0.4, 1.0);
//		randomSequence->Next();
//		g = randomSequence->GetRangeValue(0.4, 1.0);
//		randomSequence->Next();
//		b = randomSequence->GetRangeValue(0.4, 1.0);
//		randomSequence->Next();
//		actor->GetProperty()->SetDiffuseColor(r, g, b);
//		actor->GetProperty()->SetDiffuse(0.8);
//		actor->GetProperty()->SetSpecular(0.5);
//		actor->GetProperty()->SetSpecularColor(
//			colors->GetColor3d("White").GetData());
//		actor->GetProperty()->SetSpecularPower(30.0);
//		renderer->AddActor(actor);
//	}
//
//	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Initialize();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}



//vectorAcotrs sphere

//#include <vtkActor.h>
//#include <vtkInteractorStyleTrackball.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkNew.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//
//
//#include <vector>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	std::vector<vtkSmartPointer<vtkActor>> actors;
//
//	for (unsigned int i = 0; i < 10; i++)
//	{
//		vtkNew<vtkSphereSource> sphereSource;
//		sphereSource->SetCenter(i, 0.0, 0.0);
//		sphereSource->SetRadius(0.2);
//
//		vtkNew<vtkPolyDataMapper> mapper;
//		mapper->SetInputConnection(sphereSource->GetOutputPort());
//
//		vtkNew<vtkActor> actor;
//		actor->SetMapper(mapper);
//		actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
//
//		actors.push_back(actor);
//	}
//
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("VectorOfActors");
//
//	// An interactor
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the scene
//	for (unsigned int i = 0; i < actors.size(); i++)
//	{
//		renderer->AddActor(actors[i]);
//	}
//
//	renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());
//
//	// Render
//	renderWindow->Render();
//
//	vtkNew<vtkInteractorStyleTrackballCamera> style;
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	// Begin mouse interaction
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//cubeGlyph3D

//#include <vtkActor.h>
//#include <vtkCellArray.h>
//#include <vtkCubeSource.h>
//#include <vtkGlyph3D.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkPoints> points;
//	points->InsertNextPoint(0, 0, 0);
//	points->InsertNextPoint(1, 1, 0);
//	points->InsertNextPoint(2, 2, 0);
//
//	vtkNew<vtkPolyData> polydata;
//	polydata->SetPoints(points);
//
//	// Create anything you want here, we will use a cube for the demo.
//	vtkNew<vtkCubeSource> cubeSource;
//
//	vtkNew<vtkGlyph3D> glyph3D;
//	glyph3D->SetSourceConnection(cubeSource->GetOutputPort());
//	glyph3D->SetInputData(polydata);
//	glyph3D->Update();
//
//	// Visualize
//	vtkNew<vtkPolyDataMapper> mapper;
//	mapper->SetInputConnection(glyph3D->GetOutputPort());
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
//
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(actor);
//	renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());
//
//	renderWindow->SetWindowName("Glyph3D");
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}



//#include <vtkActor.h>
//#include <vtkInteractorStyleTrackball.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkNew.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkArrowSource.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//
//
//#include <vector>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	std::vector<vtkSmartPointer<vtkActor>> actors;
//
//	for (unsigned int i = 0; i < 10; i++)
//	{
//		vtkNew<vtkArrowSource> arrowSource;
//		//arrowSource->SetCenter(i, 0.0, 0.0);
//		//arrowSource->SetRadius(0.2);
//		arrowSource->set
//
//		vtkNew<vtkPolyDataMapper> mapper;
//		mapper->SetInputConnection(arrowSource->GetOutputPort());
//
//		vtkNew<vtkActor> actor;
//		actor->SetMapper(mapper);
//		actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
//
//		actors.push_back(actor);
//	}
//
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("VectorOfActors");
//
//	// An interactor
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the scene
//	for (unsigned int i = 0; i < actors.size(); i++)
//	{
//		renderer->AddActor(actors[i]);
//	}
//
//	renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());
//
//	// Render
//	renderWindow->Render();
//
//	vtkNew<vtkInteractorStyleTrackballCamera> style;
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	// Begin mouse interaction
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}



//glyph3Dmapper

//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkCellArray.h>
//#include <vtkCubeSource.h>
//#include <vtkArrowSource.h>
//#include <vtkFloatArray.h>
//#include <vtkGlyph3DMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPointData.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkUnsignedCharArray.h>
//
//int main(int, char* [])
//{
//	vtkNew<vtkPoints> points;
//	points->InsertNextPoint(0, 0, 0);
//	points->InsertNextPoint(1, 1, 1);
//	points->InsertNextPoint(2, 2, 2);
//
//	vtkNew<vtkFloatArray> scaleFactors;
//	scaleFactors->SetNumberOfComponents(3);
//	scaleFactors->SetName("Scale Factors");
//	scaleFactors->InsertNextTuple3(0.7, 1.0, 1.0);
//	scaleFactors->InsertNextTuple3(1.0, 0.7, 1.0);
//	scaleFactors->InsertNextTuple3(1.0, 1.0, 0.7);
//
//	vtkNew<vtkNamedColors> namedColors;
//
//	vtkNew<vtkUnsignedCharArray> colors;
//	colors->SetName("Colors");
//	colors->SetNumberOfComponents(3);
//	colors->InsertNextTypedTuple(namedColors->GetColor3ub("Red").GetData());
//	colors->InsertNextTypedTuple(namedColors->GetColor3ub("Green").GetData());
//	colors->InsertNextTypedTuple(namedColors->GetColor3ub("Blue").GetData());
//
//	vtkNew<vtkPolyData> polydata;
//	polydata->SetPoints(points);
//	polydata->GetPointData()->AddArray(colors);
//	polydata->GetPointData()->AddArray(scaleFactors);
//
//	// Create anything you want here, we will use a cube for the demo.
//	//vtkNew<vtkCubeSource> cubeSource;
//	vtkNew<vtkArrowSource> cubeSource;
//
//	vtkNew<vtkGlyph3DMapper> glyph3Dmapper;
//	glyph3Dmapper->SetSourceConnection(cubeSource->GetOutputPort());
//	glyph3Dmapper->SetInputData(polydata);
//	glyph3Dmapper->SetScalarModeToUsePointFieldData();
//	glyph3Dmapper->SetScaleArray("Scale Factors");
//	glyph3Dmapper->SetScaleModeToScaleByVectorComponents();
//	glyph3Dmapper->SelectColorArray("Colors");
//	glyph3Dmapper->Update();
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(glyph3Dmapper);
//
//	// Create a renderer, render window, and interactor
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("Glyph3DMapper");
//
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actor to the scene
//	renderer->AddActor(actor);
//	renderer->SetBackground(namedColors->GetColor3d("SlateGray").GetData());
//
//	// Position the camera
//	renderer->GetActiveCamera()->SetPosition(-10, 5, 0);
//	renderer->GetActiveCamera()->SetFocalPoint(1, 1, 1);
//
//	// Render and interact
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkConeSource.h>
//#include <vtkGlyph3D.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSphereSource.h>
//
//#include <array>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	// Set the background color.
//	std::array<unsigned char, 4> bkg{ {26, 51, 102, 255} };
//	colors->SetColor("Bkg", bkg.data());
//
//	// Create the rendering objects.
//	vtkNew<vtkRenderer> ren1;
//	vtkNew<vtkRenderWindow> renWin;
//	renWin->AddRenderer(ren1);
//	vtkNew<vtkRenderWindowInteractor> iren;
//	iren->SetRenderWindow(renWin);
//
//	// Create the pipeline, ball and spikes.
//	vtkNew<vtkSphereSource> sphere;
//	sphere->SetPhiResolution(7);
//	sphere->SetThetaResolution(7);
//	vtkNew<vtkPolyDataMapper> sphereMapper;
//	sphereMapper->SetInputConnection(sphere->GetOutputPort());
//	vtkNew<vtkActor> sphereActor;
//	sphereActor->SetMapper(sphereMapper);
//	vtkNew<vtkActor> sphereActor2;
//	sphereActor2->SetMapper(sphereMapper);
//
//	vtkNew<vtkConeSource> cone;
//	cone->SetResolution(5);
//	vtkNew<vtkGlyph3D> glyph;
//	glyph->SetInputConnection(sphere->GetOutputPort());
//	glyph->SetSourceConnection(cone->GetOutputPort());
//	glyph->SetVectorModeToUseNormal();
//	glyph->SetScaleModeToScaleByVector();
//	glyph->SetScaleFactor(0.25);
//	vtkNew<vtkPolyDataMapper> spikeMapper;
//	spikeMapper->SetInputConnection(glyph->GetOutputPort());
//	vtkNew<vtkActor> spikeActor;
//	spikeActor->SetMapper(spikeMapper);
//	vtkNew<vtkActor> spikeActor2;
//	spikeActor2->SetMapper(spikeMapper);
//
//	spikeActor->SetPosition(0, 0.7, 0);
//	sphereActor->SetPosition(0, 0.7, 0);
//	spikeActor2->SetPosition(0, -1.0, -10);
//	sphereActor2->SetPosition(0, -1.0, -10);
//	spikeActor2->SetScale(1.5, 1.5, 1.5);
//	sphereActor2->SetScale(1.5, 1.5, 1.5);
//
//	ren1->AddActor(sphereActor);
//	ren1->AddActor(spikeActor);
//	ren1->AddActor(sphereActor2);
//	ren1->AddActor(spikeActor2);
//	ren1->SetBackground(colors->GetColor3d("Bkg").GetData());
//	renWin->SetSize(300, 300);
//	renWin->SetWindowName("CameraBlur");
//	//   renWin->DoubleBufferOff();
//
//	// Do the first render and then zoom in a little.
//	renWin->Render();
//	ren1->GetActiveCamera()->SetFocalPoint(0, 0, 0.0);
//	ren1->GetActiveCamera()->Zoom(1.8);
//	ren1->GetActiveCamera()->SetFocalDisk(0.05);
//
//	renWin->Render();
//
//	iren->Start();
//
//	return EXIT_SUCCESS;
//}


//箭头，cellDataToPoint，单元数据转成节点数据，选择标量显示

//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkPolyData.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkLookupTable.h>
//#include <algorithm>
//#include <array>
//#include <string>
//#include <vtkScalarBarActor.h>
//#include <vtkAxesActor.h>
//#include <vtkOrientationMarkerWidget.h>
//#include <vtkDataSetAttributes.h>
//#include <vtkPolyData.h>
//#include <vtkPointData.h>
//#include <vtkCellData.h>
//#include <vtkArrowSource.h>
//#include <vtkGlyph3D.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkGlyphSource2D.h>
//#include <vtkCellDataToPointData.h>
//#include <vtkCellData.h>
//#include <vtkMaskPoints.h>
//#include <vtkUnstructuredGrid.h>
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//
//	vtkNew<vtkNamedColors> colors;
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//
//
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkInteractorStyleTrackballCamera> style;
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//
//	interactor->SetInteractorStyle(style);
//	interactor->SetRenderWindow(renderWindow);
//
//	renderer->SetBackground(colors->GetColor3d("LightSlateGray").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//vtkNew<vtkXMLUnstructuredGridReader>reader;
//	reader->ReadAllScalarsOn();//获取所有的标量数据
//	reader->ReadAllVectorsOn();
//	reader->ReadAllNormalsOn();
//	reader->ReadAllTensorsOn();
//	reader->ReadAllColorScalarsOn();
//	reader->ReadAllTCoordsOn();
//	reader->ReadAllFieldsOn();
//	
//
//	//reader->SetFileName("step22.vtk");
//	reader->SetFileName("1.vtk");
//	//reader->GetOutput()->Register(reader);
//	reader->Update();
//
//	vtkNew<vtkCellDataToPointData>cellToPoint;
//	cellToPoint->SetInputData(reader->GetOutput());
//	cellToPoint->Update();
//
//	vtkNew<vtkMaskPoints>mask;
//	mask->SetInputData(cellToPoint->GetOutput());
//	mask->SetMaximumNumberOfPoints(9000);
//	mask->RandomModeOn();
//	mask->Update();
//
//	int cellNum = cellToPoint->GetOutput()->GetNumberOfCells();
//	cout << "cells of cellToPoint: " << cellNum << endl;
//	int pointNum = cellToPoint->GetOutput()->GetNumberOfPoints();
//	cout << "points of cellToPoint: " << pointNum<< endl;
//
//	int cellNum2 = reader->GetOutput()->GetNumberOfCells();
//	cout << "cells of reader: " << cellNum2 << endl;
//	int pointNum2 = reader->GetOutput()->GetNumberOfPoints();
//	cout << "points of reader: " << pointNum2 << endl;
//
//	int nNumScalar = reader->GetNumberOfScalarsInFile();//获取标量类型数
//	cout << "Number of scalars is: " << nNumScalar << endl;
//
//	vtkSmartPointer<vtkUnstructuredGrid> UnstructuredGrid =
//		vtkSmartPointer<vtkUnstructuredGrid>::New();
//	UnstructuredGrid = reader->GetOutput();
//	auto unstructuredGrid = reader->GetOutput();
//	for (int i = 0; i < nNumScalar; i++)
//	{
//		cout << reader->GetScalarsNameInFile(i) << endl;
//	}
//
//	int nNumVector = reader->GetNumberOfVectorsInFile();
//	cout << "Number of Vector is: " << nNumVector << endl;
//	
//	for (int i = 0; i < nNumVector; i++)
//	{
//		cout << reader->GetVectorsNameInFile(i) << endl;
//	}
//
//	//vtkNew<vtkArrowSource> arrow;
//	vtkNew<vtkGlyphSource2D>arrow;//二维箭头
//	arrow->SetGlyphTypeToArrow();
//	arrow->FilledOff();
//	
//	vtkNew<vtkGlyph3D> glyphs;
//	glyphs->SetSourceConnection(arrow->GetOutputPort());
//	glyphs->SetInputConnection(mask->GetOutputPort());
//	//glyphs->SetInputConnection(cellToPoint->GetOutputPort());
//	glyphs->ScalingOn();
//	glyphs->SetScaleModeToScaleByVector();
//	glyphs->SetScaleFactor(0.005);
//	glyphs->OrientOn();
//	glyphs->ClampingOff();
//	glyphs->SetVectorModeToUseVector();
//	glyphs->SetIndexModeToOff();
//	//glyphs->SetColorModeToColorByVector();
//	
//	int scalarIndex = 4;
//
//	reader->GetOutput()->GetPointData()->SetActiveScalars(reader->GetScalarsNameInFile(scalarIndex));//设置标量名称，即渲染哪个标量
//	cellToPoint->GetOutput()->GetPointData()->SetActiveScalars(reader->GetScalarsNameInFile(scalarIndex));//设置标量名称，即渲染哪个标量
//	//reader->GetOutput()->GetCellData()->SetActiveScalars(reader->GetScalarsNameInFile(scalarIndex));
//	reader->Update();
//
//	vtkNew<vtkLookupTable> lut1;
//	//lut1->SetHueRange(0.5, 0.833);// 设定HSV颜色范围，色调H取值范围为0°～360°，从红色开始按逆时针方向计算，红色为0°/0.0，绿色为120°/0.34,蓝色为240°/0.67
//	lut1->SetHueRange(.667, 0);
//
//	vtkNew<vtkPolyDataMapper> glyphMapper;
//	glyphMapper->SetInputConnection(glyphs->GetOutputPort());
//	//glyphMapper->SetScalarRange(cellToPoint->GetOutput()->GetScalarRange());
//	glyphMapper->SetLookupTable(lut1);
//	glyphMapper->SetColorModeToMapScalars();
//
//	// Visualize
//	vtkNew<vtkDataSetMapper> mapper;
//	//mapper->SetInputData(unstructuredGrid);
//	mapper->SetInputData(cellToPoint->GetOutput());
//	//mapper->SetScalarRange(unstructuredGrid->GetScalarRange());
//	mapper->SetScalarRange(cellToPoint->GetOutput()->GetScalarRange());
//
//	mapper->SetLookupTable(lut1);
//	mapper->SetColorModeToMapScalars();
//	//mapper->SetScalarModeToUseCellData();//使用网格数据进行映射
//	mapper->SetScalarModeToDefault();
//
//	cout << "unstructuredGrid's scalar range is: " << endl;
//	cout << unstructuredGrid->GetScalarRange()[0] << endl;
//	cout << unstructuredGrid->GetScalarRange()[1] << endl;
//
//	vtkNew<vtkScalarBarActor> scalarbar;
//	scalarbar->SetLookupTable(glyphMapper->GetLookupTable());
//	//scalarbar->SetTitle(curvaturesfilter->GetOutput()->GetPointData()->GetScalars()->GetName());
//	scalarbar->SetNumberOfLabels(8);
//	scalarbar->SetMaximumWidthInPixels(100);
//	//scalarbar->SetTitle("Velocity(m/s)");
//	scalarbar->SetTitle(reader->GetScalarsNameInFile(scalarIndex));
//	//scalarbar->SetTitle(reader->GetScalarsName());
//
//	renderer->AddActor2D(scalarbar);
//
//	vtkNew<vtkAxesActor> axes;
//
//	vtkNew<vtkOrientationMarkerWidget> widget;
//	double rgba[4]{ 0.0, 0.0, 0.0, 0.0 };
//	colors->GetColor("Carrot", rgba);
//	widget->SetOutlineColor(rgba[0], rgba[1], rgba[2]);
//	widget->SetOrientationMarker(axes);
//	widget->SetInteractor(interactor);
//	widget->SetViewport(0.0, 0.0, 0.2, 0.2);
//	widget->SetEnabled(1);
//	//固定坐标轴位置
//	//widget->InteractiveOn();
//	widget->InteractiveOff();
//	widget->KeyPressActivationOff();
//	widget->Modified();
//
//
//	vtkNew<vtkActor> glyphActor;
//	glyphActor->SetMapper(glyphMapper);
//	glyphActor->GetProperty()->EdgeVisibilityOff();//不显示网格
//	glyphActor->GetProperty()->SetOpacity(1);
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(mapper);
//	//actor->GetProperty()->EdgeVisibilityOn();//显示网格
//	actor->GetProperty()->SetOpacity(1);//透明度
//	//actor->GetProperty()->SetRepresentationToWireframe();//显示外网格
//
//	//renderer->AddActor(glyphActor);
//	renderer->AddActor(actor);
//
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}



//generate point normals using local tangent planes 指向球体外法向的球体轮廓箭头

//#include <vtkArrowSource.h>
//#include <vtkCamera.h>
//#include <vtkGlyph3D.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPCANormalEstimation.h>
//#include <vtkPointSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSphereSource.h>
//
//namespace {
//	void MakeGlyphs(vtkPolyData* src, double size, vtkGlyph3D* glyph);
//}
//
//int main(int, char* [])
//{
//	double radius = 1.0;
//	vtkNew<vtkPointSource> points;
//	points->SetNumberOfPoints(600);
//	points->SetRadius(radius);
//	points->SetCenter(0.0, 0.0, 0.0);
//	points->SetDistributionToShell();
//
//	int sampleSize = 10;
//	vtkNew<vtkPCANormalEstimation> normals;
//	normals->SetInputConnection(points->GetOutputPort());
//	normals->SetSampleSize(sampleSize);
//	normals->SetNormalOrientationToGraphTraversal();
//	normals->Update();
//
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkGlyph3D> glyph3D;
//	MakeGlyphs(normals->GetOutput(), radius * 0.3, glyph3D.GetPointer());
//
//	vtkNew<vtkPolyDataMapper> glyph3DMapper;
//	glyph3DMapper->SetInputConnection(glyph3D->GetOutputPort());
//
//	vtkNew<vtkActor> glyph3DActor;
//	glyph3DActor->SetMapper(glyph3DMapper);
//	glyph3DActor->GetProperty()->SetDiffuseColor(
//		colors->GetColor3d("Banana").GetData());
//
//	vtkNew<vtkSphereSource> sphere;
//	sphere->SetRadius(1.0);
//	sphere->SetThetaResolution(41);
//	sphere->SetPhiResolution(21);
//
//	vtkNew<vtkPolyDataMapper> sphereMapper;
//	sphereMapper->SetInputConnection(sphere->GetOutputPort());
//
//	vtkNew<vtkActor> sphereActor;
//	sphereActor->SetMapper(sphereMapper);
//	sphereActor->GetProperty()->SetDiffuseColor(
//		colors->GetColor3d("Tomato").GetData());
//
//	// Create graphics stuff
//	//
//	vtkNew<vtkRenderer> renderer;
//	renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetSize(640, 480);
//	renderWindow->SetWindowName("NormalEstimation");
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the renderer, set the background and size
//	//
//	renderer->AddActor(glyph3DActor);
//	renderer->AddActor(sphereActor);
//
//	// Generate an interesting view
//	//
//	renderer->ResetCamera();
//	renderer->GetActiveCamera()->Azimuth(120);
//	renderer->GetActiveCamera()->Elevation(30);
//	renderer->GetActiveCamera()->Dolly(1.0);
//	renderer->ResetCameraClippingRange();
//
//	renderWindow->Render();
//	interactor->Initialize();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}
//namespace {
//	void MakeGlyphs(vtkPolyData* src, double size, vtkGlyph3D* glyph)
//	{
//		// Source for the glyph filter
//		vtkNew<vtkArrowSource> arrow;
//		arrow->SetTipResolution(16);
//		arrow->SetTipLength(0.3);
//		arrow->SetTipRadius(0.1);
//
//		glyph->SetSourceConnection(arrow->GetOutputPort());
//		glyph->SetInputData(src);
//		glyph->SetVectorModeToUseNormal();
//		glyph->SetScaleModeToScaleByVector();
//		glyph->SetScaleFactor(size);
//		glyph->OrientOn();
//		glyph->Update();
//	}
//} // namespace



//艺术样条放样（有意思）

//#include <vtkActor.h>
//#include <vtkCardinalSpline.h>
//#include <vtkCellArray.h>
//#include <vtkGlyph3D.h>
//#include <vtkMath.h>
//#include <vtkMinimalStandardRandomSequence.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPoints.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkSplineFilter.h>
//#include <vtkXMLPolyDataReader.h>
//
//int main(int argc, char* argv[])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	auto polyData = vtkSmartPointer<vtkPolyData>::New();
//	if (argc > 1)
//	{
//		vtkNew<vtkXMLPolyDataReader> reader;
//		reader->SetFileName(argv[1]);
//		reader->Update();
//		polyData = reader->GetOutput();
//	}
//	else
//	{
//		unsigned int numberOfPoints = 100;
//		vtkNew<vtkPoints> points;
//		vtkNew<vtkMinimalStandardRandomSequence> randomSequence;
//		randomSequence->SetSeed(8775070);
//		for (unsigned int i = 0; i < numberOfPoints; ++i)
//		{
//			double x, y, z;
//			// random position and radius
//			x = randomSequence->GetRangeValue(-1.0, 1.0);
//			randomSequence->Next();
//			y = randomSequence->GetRangeValue(-1.0, 1.0);
//			randomSequence->Next();
//			z = randomSequence->GetRangeValue(-1.0, 1.0);
//			randomSequence->Next();
//			points->InsertNextPoint(x, y, z);
//		}
//		vtkNew<vtkCellArray> lines;
//		lines->InsertNextCell(numberOfPoints);
//		for (unsigned int i = 0; i < numberOfPoints; ++i)
//		{
//			lines->InsertCellPoint(i);
//		}
//		polyData->SetPoints(points);
//		polyData->SetLines(lines);
//	}
//
//	vtkNew<vtkCardinalSpline> spline;
//	spline->SetLeftConstraint(2);
//	spline->SetLeftValue(0.0);
//	spline->SetRightConstraint(2);
//	spline->SetRightValue(0.0);
//
//	vtkNew<vtkSplineFilter> splineFilter;
//	splineFilter->SetInputData(polyData);
//	splineFilter->SetNumberOfSubdivisions(polyData->GetNumberOfPoints() * 10);
//	splineFilter->SetSpline(spline);
//
//	vtkNew<vtkPolyDataMapper> splineMapper;
//	splineMapper->SetInputConnection(splineFilter->GetOutputPort());
//
//	vtkNew<vtkActor> splineActor;
//	splineActor->SetMapper(splineMapper);
//
//	vtkNew<vtkSphereSource> originalNodes;
//	originalNodes->SetRadius(.04);
//	originalNodes->SetPhiResolution(10);
//	originalNodes->SetThetaResolution(10);
//
//	vtkNew<vtkGlyph3D> glyphOriginal;
//	glyphOriginal->SetInputData(polyData);
//	glyphOriginal->SetSourceConnection(originalNodes->GetOutputPort());
//
//	vtkNew<vtkSphereSource> newNodes;
//	newNodes->SetRadius(.02);
//	newNodes->SetPhiResolution(10);
//	newNodes->SetThetaResolution(10);
//
//	vtkNew<vtkGlyph3D> glyphNew;
//	glyphNew->SetInputConnection(splineFilter->GetOutputPort());
//	glyphNew->SetSourceConnection(newNodes->GetOutputPort());
//
//	vtkNew<vtkPolyDataMapper> originalMapper;
//	originalMapper->SetInputConnection(glyphOriginal->GetOutputPort());
//
//	vtkNew<vtkActor> originalActor;
//	originalActor->SetMapper(originalMapper);
//	originalActor->GetProperty()->SetColor(
//		colors->GetColor3d("Banana").GetData());
//	originalActor->GetProperty()->SetOpacity(.6);
//
//	vtkNew<vtkPolyDataMapper> newMapper;
//	newMapper->SetInputConnection(glyphNew->GetOutputPort());
//
//	vtkNew<vtkActor> newActor;
//	newActor->SetMapper(newMapper);
//	newActor->GetProperty()->SetColor(colors->GetColor3d("Tomato").GetData());
//
//	// A renderer and render window
//	vtkNew<vtkRenderer> renderer;
//	renderer->SetBackground(colors->GetColor3d("SteelBlue").GetData());
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("ResamplePolyLine");
//
//	// An interactor
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actors to the scene
//	renderer->AddActor(originalActor);
//	renderer->AddActor(newActor);
//	renderer->AddActor(splineActor);
//
//	renderWindow->Render();
//
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}



//VisualizeDirectedGraph 首尾相接形成三角形的箭头

//#include <vtkActor.h>
//#include <vtkGlyph3D.h>
//#include <vtkGlyphSource2D.h>
//#include <vtkGraphLayout.h>
//#include <vtkGraphLayoutView.h>
//#include <vtkGraphToPolyData.h>
//#include <vtkMutableDirectedGraph.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSimple2DLayoutStrategy.h>
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkMutableDirectedGraph> g;
//
//	vtkIdType v1 = g->AddVertex();
//	vtkIdType v2 = g->AddVertex();
//	vtkIdType v3 = g->AddVertex();
//
//	g->AddEdge(v1, v2);
//	g->AddEdge(v2, v3);
//	g->AddEdge(v3, v1);
//
//	// Do layout manually before handing graph to the view.
//	// This allows us to know the positions of edge arrows.
//	vtkNew<vtkGraphLayoutView> graphLayoutView;
//
//	vtkNew<vtkGraphLayout> layout;
//	vtkNew<vtkSimple2DLayoutStrategy> strategy;
//	layout->SetInputData(g);
//	layout->SetLayoutStrategy(strategy);
//
//	// Tell the view to use the vertex layout we provide
//	graphLayoutView->SetLayoutStrategyToPassThrough();
//	// The arrows will be positioned on a straight line between two
//	// vertices so tell the view not to draw arcs for parallel edges
//	graphLayoutView->SetEdgeLayoutStrategyToPassThrough();
//
//	// Add the graph to the view. This will render vertices and edges,
//	// but not edge arrows.
//	graphLayoutView->AddRepresentationFromInputConnection(
//		layout->GetOutputPort());
//
//	// Manually create an actor containing the glyphed arrows.
//	vtkNew<vtkGraphToPolyData> graphToPoly;
//	graphToPoly->SetInputConnection(layout->GetOutputPort());
//	graphToPoly->EdgeGlyphOutputOn();
//
//	// Set the position (0: edge start, 1: edge end) where
//	// the edge arrows should go.
//	graphToPoly->SetEdgeGlyphPosition(0.98);
//
//	// Make a simple edge arrow for glyphing.
//	vtkNew<vtkGlyphSource2D> arrowSource;
//	arrowSource->SetGlyphTypeToEdgeArrow();
//	arrowSource->SetScale(0.1);
//	arrowSource->Update();
//
//	// Use Glyph3D to repeat the glyph on all edges.
//	vtkNew<vtkGlyph3D> arrowGlyph;
//	arrowGlyph->SetInputConnection(0, graphToPoly->GetOutputPort(1));
//	arrowGlyph->SetInputConnection(1, arrowSource->GetOutputPort());
//
//	// Add the edge arrow actor to the view.
//	vtkNew<vtkPolyDataMapper> arrowMapper;
//	arrowMapper->SetInputConnection(arrowGlyph->GetOutputPort());
//	vtkNew<vtkActor> arrowActor;
//	arrowActor->SetMapper(arrowMapper);
//	graphLayoutView->GetRenderer()->AddActor(arrowActor);
//
//	graphLayoutView->GetRenderer()->SetBackground(
//		colors->GetColor3d("SaddleBrown").GetData());
//	graphLayoutView->GetRenderer()->SetBackground2(
//		colors->GetColor3d("Wheat").GetData());
//	graphLayoutView->GetRenderWindow()->SetWindowName("VisualizeDirectedGraph");
//	graphLayoutView->ResetCamera();
//	graphLayoutView->Render();
//	graphLayoutView->GetInteractor()->Start();
//
//	return EXIT_SUCCESS;
//}



//vtu转换为vtk文件格式

//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <vtkXMLPolyDataReader.h>
//#include <vtkPLYWriter.h>
//#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkUnstructuredGridWriter.h>
//
//int main(int argc, char* argv[])
//{
//	for (int i = 9; i < 51; i++)
//	{
//		//std::string inputFileName = "vtuReaderTest.vtu";
//		int fileNum = i;
//		std::string number = to_string(i + 1);
//		std::string inputFileName = "djc.t0" + number;
//		inputFileName += ".p000.vtu";
//
//		std::string outputFileName = "step" + number;
//		outputFileName += ".vtk";
//		cout << outputFileName << endl;
//
//		vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
//		reader->SetFileName(inputFileName.c_str());
//		reader->Update();
//
//		vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
//		writer->SetFileName(outputFileName.c_str());
//		writer->SetInputConnection(reader->GetOutputPort());
//		writer->Update();
//		cout << outputFileName << " is finished." << endl;
//	}
//	return 0;
//}


//法线 二维四边形单元

//#include <cmath>
//#include <iostream>
//#include <vector>
//
//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataNormals.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkCellData.h>
//#include <vtkPointData.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkConeSource.h>
//#include <vtkTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkGlyph3D.h>
//#include <vtkGlyphSource2D.h>
//#include <vtkCellCenters.h>
//
//#include <algorithm>
//#include <array>
//#include <string>
//#include <fstream>
//// #include <iostream>
//
//
//class Vector3
//{
//public:
//	Vector3();
//	Vector3(double X_, double Y_, double Z_)
//	{
//		m_x = X_;
//		m_y = Y_;
//		m_z = Z_;
//	}
//
//	const double x() { return m_x; }
//	const double y() { return m_y; }
//	const double z() { return m_z; }
//
//	Vector3 cross(Vector3& vec)
//	{
//		double X_ = this->y() * vec.z() - this->z() * vec.y();
//		double Y_ = this->z() * vec.x() - this->x() * vec.z();
//		double Z_ = this->x() * vec.y() - this->y() * vec.x();
//
//		return Vector3(X_, Y_, Z_);
//	}
//
//	void normalized()
//	{
//		double w = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
//		m_x /= w;
//		m_y /= w;
//		m_z /= w;
//	}
//
//private:
//	double m_x;
//	double m_y;
//	double m_z;
//};
//
////给定三个点p1 p2 p3，返回由这三个点确定平面的法向n，已归一化
//Vector3 getNormal(vector<double> p1, vector<double> p2, vector<double> p3)
//{
//	Vector3 v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
//	Vector3 v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
//
//	Vector3 n = v1.cross(v2);
//	n.normalized();
//
//	return n;
//}
//
//vector<double> normalized(vector<double>& vector)
//{
//	double L = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
//	vector[0] /= L;
//	vector[1] /= L;
//	vector[2] /= L;
//	return vector;
//}
//
///*zNormal:该单元法向
//* kSet:用户设置的x或y方向
//* point:该单元中心坐标
//* 返回x或y方向渗透率的投影方向
//*/
//vector<double>  projectionNormal(vector<double> zNormal, const vector<double> kSet, const double* point)
//{
//	//分别计算xNor（kSet[0]）和yNor（kSet[1]）
//	double Q[3];
//	for (int i = 0; i < 3; i++)
//	{
//		Q[i] = point[i] + kSet[i];
//	}
//
//	//temp: k在n方向上的投影数值大小
//	double temp = (kSet[0] * zNormal[0] + kSet[1] * zNormal[1] + kSet[2] * zNormal[2]) /
//		sqrt(zNormal[0] * zNormal[0] + zNormal[1] * zNormal[1] + zNormal[2] * zNormal[2]);
//
//	//单位化zNormal
//	double n[3];
//	n[0] = normalized(zNormal)[0];
//	n[1] = normalized(zNormal)[1];
//	n[2] = normalized(zNormal)[2];
//
//	double j[3];
//	for (int i = 0; i < 3; i++)
//	{
//		j[i] = temp * n[i];
//	}
//
//	//R: 投影后的点R
//	double R[3];
//	for (int i = 0; i < 3; i++)
//	{
//		R[i] = Q[i] - j[i];
//	}
//
//	//r: 投影后向量方向
//	double r[3];
//	for (int i = 0; i < 3; i++)
//	{
//		r[i] = R[i] - point[i];
//	}
//
//	//结果写入projNor
//	vector<double> projNor(3);
//	for (int i = 0; i < 3; i++)
//	{
//		projNor[i] = r[i];
//	}
//
//	return projNor;
//}
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	//renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->SetBackground(colors->GetColor3d("Black").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//reader->SetFileName("Test01.vtk");
//	//reader->SetFileName("step12.vtk");
//	//reader->SetFileName("blade.vtk");
//	//reader->SetFileName("plane.vtk");
//	//reader->SetFileName("S.vtk");
//	reader->SetFileName("qu.vtk");
//	//reader->SetFileName("bladeGmsh3_vtk.vtk");
//	reader->Update();
//
//	auto unstructuredGrid = reader->GetOutput();
//
//	// 叉乘得到的法向
//	vector<double> calNor;
//
//	// 单元中心坐标
//	vector<double> centrePoint;
//	
//	for (int i = 0; i < unstructuredGrid->GetNumberOfCells(); i++)
//	{
//		vector<double> p1, p2, p3, p4;
//
//		vtkPoints* point0 = unstructuredGrid->GetCell(i)->GetPoints();
//		double* coors1 = point0->GetPoint(0);
//		p1.push_back(*coors1);
//		p1.push_back(*(coors1 + 1));
//		p1.push_back(*(coors1 + 2));
//
//		double* coors2 = point0->GetPoint(1);
//		p2.push_back(*coors2);
//		p2.push_back(*(coors2 + 1));
//		p2.push_back(*(coors2 + 2));
//
//		double* coors3 = point0->GetPoint(2);
//		p3.push_back(*coors3);
//		p3.push_back(*(coors3 + 1));
//		p3.push_back(*(coors3 + 2));
//
//		double* coors4 = point0->GetPoint(3);
//		p4.push_back(*coors4);
//		p4.push_back(*(coors4 + 1));
//		p4.push_back(*(coors4 + 2));
//
//		Vector3 cal_nor = getNormal(p1, p2, p3);
//
//		// 把p1 p2 p3平面的法向按顺序存进calNor
//		calNor.push_back(cal_nor.x());
//		calNor.push_back(cal_nor.y());
//		calNor.push_back(cal_nor.z());
//
//		double x = (p1[0] + p2[0] + p3[0] + p4[0]) / 4;
//		double y = (p1[1] + p2[1] + p3[1] + p4[1]) / 4;
//		double z = (p1[2] + p2[2] + p3[2] + p4[2]) / 4;
//
//
//		// 把i单元的中心点按顺序存进centrePoint
//		centrePoint.push_back(x);
//		centrePoint.push_back(y);
//		centrePoint.push_back(z);
//	}
//
//	ofstream ofs2;
//	ofs2.open("calculate_Normal.txt", ios::out);
//	ofs2 << "For points ---------------------: \n";
//
//	for (int i = 0; i < calNor.size(); i+=3)
//	{
//		ofs2 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	// unstructuredGrid 转换成 polydata
//	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
//		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//	surfaceFilter->SetInputData(unstructuredGrid);
//	surfaceFilter->Update();
//
//	vtkPolyData* polydata = surfaceFilter->GetOutput();
//
//	/*------------------------------数据可视化------------------------------*/
//	// Visualize
//	vtkNew<vtkDataSetMapper> dataSetMapper;
//	dataSetMapper->SetInputData(unstructuredGrid);
//	dataSetMapper->ScalarVisibilityOff();
//	
//	vtkNew<vtkProperty> backProp;
//	backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
//	backProp->SetSpecular(.6);
//	backProp->SetSpecularPower(30);
//	
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(dataSetMapper);
//	actor->SetBackfaceProperty(backProp);
//	actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
//	actor->GetProperty()->SetSpecular(.3);
//	actor->GetProperty()->SetSpecularPower(30);
//	actor->GetProperty()->EdgeVisibilityOn();
//	renderer->AddActor(actor);
//	//
//	/*renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();*/
//	/*-----------------------------------------------------------------*/
//
//	vtkSmartPointer<vtkPolyDataMapper> mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(surfaceFilter->GetOutputPort());
//
//	// 测试unstructuredgrid转为surface后的坐标顺序 ： 以单元编号开始输出对应的点坐标
//	//for (int i = 0; i < 200; i++)
//	//{
//	//	double* point = surfaceFilter->GetOutput()->GetPoint(i);
//	//	cout << i << " " << point[0] << " " << point[1] << " " << point[2] << endl;
//	//}
//
//	vtkSmartPointer<vtkActor> surfaceActor =
//		vtkSmartPointer<vtkActor>::New();
//	surfaceActor->SetMapper(mapper);
//	vtkSmartPointer<vtkPolyDataNormals> pdNormals =
//		vtkSmartPointer<vtkPolyDataNormals>::New();
//	pdNormals->SetInputConnection(surfaceFilter->GetOutputPort());
//	pdNormals->ComputeCellNormalsOn();
//	pdNormals->Update();
//
//	vtkPointData* ptData = pdNormals->GetOutput()->GetPointData();
//
//	ofstream ofs;
//	vector<double*> zNormal;
//	vector<double> kx;
//	vector<double> ky;
//	vector<double> kz;
//	vector<double> cellKz;
//
//	//两个方向kx和ky
//	vector<double> kxSet{ 1.0, 0.0, 0.0 };
//	vector<double> kySet{ 0.0, 1.0, 0.0 };
//
//	if (ptData)
//	{
//		vtkDataArray* ptNormals = pdNormals->GetOutput()->GetPointData()->GetNormals();
//		if (ptNormals)
//		{
//			//ofs创建normal.txt，写入文件头
//			ofs.open("normal.txt", ios::out);
//			ofs << "For points in every cell: \n";
//			ofs << ptNormals->GetNumberOfTuples() << endl;
//
//			cout << "For points in every cell: \n";
//			cout << ptNormals->GetNumberOfTuples() << endl;
//			for (int i = 0; i < ptNormals->GetNumberOfTuples(); ++i)
//			{
//				double value[3];	//索引为i的法线方向
//				ptNormals->GetTuple(i, value);
//				
//				// 输出点上法线方向
//				ofs << i << " " << value[0] << " " << value[1] << " " << value[2] << endl;
//
//				zNormal.push_back(value);	//装入zNormal
//
//				// 是否换成kx、ky、kz ???  是
//				kz.push_back(*value);
//				kz.push_back(*(value + 1));
//				kz.push_back(*(value + 2));
//
//				//	测试zNormal是否正确：正确
//				std::cout << "ZNORMAL COUT::  " << (double)zNormal[i][0] << " " << (double)zNormal[i][1] << " " << (double)zNormal[i][2] << std::endl;
//				cout << "KZ COUT:  " << kz[(0 + 3 * i)] << " " << kz[(1 + 3 * i)] << " " << kz[(2 + 3 * i)] << endl;
//				cout << endl;
//			}
//
//			//测试循环外Kz数值是否正确：正确
//			for (int k = 0; k < kz.size(); k += 3)
//			{
//				cout << "AFTER FOR KZ COUT:  " << kz[k] << " " << kz[k + 1] << " " << kz[k + 2] << endl;
//			}
//		}
//	}
//	//单元上法线
//
//	//ofs << endl;
//	//ofs << "For cells: \n";
//	//cout << "For cells: \n";
//	//if (pdNormals->GetOutput()->GetCellData() && pdNormals->GetOutput()->GetCellData()->GetNormals())
//	//{
//	//	vtkDataArray* cellNormals = pdNormals->GetOutput()->GetCellData()->GetNormals();
//	//
//	//	ofs.open("normal.txt", ios::out);
//	//	ofs << "For points in every cell: \n";
//	//	ofs << cellNormals->GetNumberOfTuples() << endl;
//	//
//	//	cout << "For points in every cell: \n";
//	//	cout << cellNormals->GetNumberOfTuples() << endl;
//	//	for (int i = 0; i < cellNormals->GetNumberOfTuples(); ++i)
//	//	{
//	//		double value[3];	//索引为i的法线方向
//	//		cellNormals->GetTuple(i, value);
//	//		printf("Value: (%d, %lf, %lf, %lf)\n", i, value[0], value[1], value[2]);
//	//		ofs << i << " " << value[0] << " " << value[1] << " " << value[2] << endl;
//	//
//	//		zNormal.push_back(value);	//装入zNormal
//	//
//	//		// 是否换成kx、ky、kz ???  是
//	//		kz.push_back(*value);
//	//		kz.push_back(*(value + 1));
//	//		kz.push_back(*(value + 2));
//	//
//	//		//	测试zNormal是否正确：正确
//	//		printf("zNormal: (%d, %lf, %lf, %lf)\n", i, zNormal[i][0], zNormal[i][1], zNormal[i][2]);
//	//		std::cout << "ZNORMAL COUT::  " << (double)zNormal[i][0] << " " << (double)zNormal[i][1] << " " << (double)zNormal[i][2] << std::endl;
//	//		cout << "KZ COUT:  " << kz[(0 + 3 * i)] << " " << kz[(1 + 3 * i)] << " " << kz[(2 + 3 * i)] << endl;
//	//		cout << endl;
//	//	}
//
//	ofstream ofs1;
//	ofs1.open("normalAndPoint.txt", ios::out);
//	ofs1 << "For points ---------------------: \n";
//
//	// 投影
//	for (int i = 0; i < calNor.size(); i+=3)
//	{
//		// double* point = pdNormals->GetOutput()->GetPoint(i);
//		double point[3] = { centrePoint[i], centrePoint[i + 1], centrePoint[i + 2] };
//		vector<double> n{ calNor[i], calNor[i + 1], calNor[i + 2] };
//
//		ofs1 << point[0] << " " << point[1] << " " << point[2] << endl;
//
//		vector<double> xProjNor(3);	//x方向上的投影
//
//		xProjNor = projectionNormal(n, kxSet, point);
//		// 计算出的kxSet放入kx中
//		kx.push_back(xProjNor[0]);
//		kx.push_back(xProjNor[1]);
//		kx.push_back(xProjNor[2]);
//
//		vector<double> yProjNor(3);	// y方向上投影
//
//		yProjNor = projectionNormal(n, kySet, point);
//		// 计算出的kySet放入ky中
//		ky.push_back(yProjNor[0]);
//		ky.push_back(yProjNor[1]);
//		ky.push_back(yProjNor[2]);
//	}
//
//	ofs1 << "For Normals  kx  ---------------------: \n";
//	for (int i = 0; i < kx.size(); i+=3)
//	{
//		ofs1 << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  ky  ---------------------: \n";
//	for (int i = 0; i < ky.size(); i+=3)
//	{
//		ofs1 << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  kz  ---------------------: \n";
//	for (int i = 0; i < calNor.size(); i+=3)
//	{
//		ofs1 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	// 单元法向量获取
//	if (pdNormals->GetOutput()->GetCellData() && pdNormals->GetOutput()->GetCellData()->GetNormals())
//	{
//		vtkDataArray* cellNormals = pdNormals->GetOutput()->GetCellData()->GetNormals();
//
//		ofs << "For cell number: \n";
//		ofs << cellNormals->GetNumberOfTuples() << endl;
//
//		cout << "For points in every cell: \n";
//		cout << cellNormals->GetNumberOfTuples() << endl;
//		for (int i = 0; i < cellNormals->GetNumberOfTuples(); ++i)
//		{
//			double value[3];	//索引为i的法线方向
//			cellNormals->GetTuple(i, value);
//
//			cellKz.push_back(*value);
//			cellKz.push_back(*(value + 1));
//			cellKz.push_back(*(value + 2));
//		}
//	}
//
//	ofs1 << "For Normals  CellKz  ---------------------: \n";
//	for (int i = 0; i < cellKz.size(); i+=3)
//	{
//		ofs1 << cellKz[i] << " " << cellKz[i + 1] << " " << cellKz[i + 2] << endl;
//	}
//
//	for (int i = 0; i < kx.size(); i += 3)
//	{
//		ofs << i / 3 << " " << "kx: " << " " << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	for (int i = 0; i < ky.size(); i += 3)
//	{
//		ofs << i / 3 << " " << "ky: " << " " << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//
//
//	vtkNew<vtkGlyphSource2D>arrow;//二维箭头
//	arrow->SetGlyphTypeToArrow();
//	arrow->FilledOff();
//
//	vtkSmartPointer<vtkTransform> transform =
//		vtkSmartPointer<vtkTransform>::New();
//	transform->RotateY(180); // make vertex outside
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformF =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	//transformF->SetInputConnection(cone->GetOutputPort());
//	transformF->SetInputConnection(arrow->GetOutputPort());
//	transformF->SetTransform(transform);
//
//	vtkSmartPointer<vtkGlyph3D> glyph =
//		vtkSmartPointer<vtkGlyph3D>::New();
//	glyph->SetInputConnection(pdNormals->GetOutputPort());
//	//glyph->SetSourceConnection(transformF->GetOutputPort()); // source => transform => graph3D
//	glyph->SetSourceConnection(arrow->GetOutputPort()); // source => transform => graph3D
//	glyph->SetVectorModeToUseNormal();
//	glyph->SetScaleModeToScaleByVector();
//	glyph->SetScaleFactor(16);
//
//	vtkSmartPointer<vtkPolyDataMapper> spikeMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	spikeMapper->SetInputConnection(glyph->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> spikeActor = vtkSmartPointer<vtkActor>::New();
//	spikeActor->SetMapper(spikeMapper);
//	spikeActor->GetProperty()->SetColor(0.0, 0.79, 0.34);
//
//	// Add the actors to the renderer, set the background and size
//	//renderer->AddActor(surfaceActor);
//	renderer->AddActor(spikeActor);
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}



//设置方向

//#include <cmath>
//#include <iostream>
//#include <vector>
//
//class HsPoint
//{
//public:
//	HsPoint(double X_, double Y_, double Z_)
//	{
//		x = X_;
//		y = Y_;
//		z = Z_;
//	}
//	~HsPoint() {};
//
//private:
//	double x;
//	double y;
//	double z;
//};
//
//double* normalized(double* vector)
//{
//	double L = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
//	vector[0] /= L;
//	vector[1] /= L;
//	vector[2] /= L;
//	return vector;
//}
//
///*zNormal:单元法向
//* kSet:用户设置的x和y方向
//* point:单元中心坐标
//* 返回3 x 3的矩阵
//*/
//vector<vector<double>>  projectionNormal(double* zNormal, vector<vector<double>> kSet, double* point)
//{
//	vector<vector<double>> projNor(3, vector<double>(3));
//	
//	//写入z方向
//	for (int i = 0; i < 3; i++)
//	{
//		projNor[2][i] = zNormal[i];
//	}
//
//	//分别计算xNor（kSet[0]）和yNor（kSet[1]）
//	for (int i = 0; i < 2; i++)
//	{
//		double Q[3];
//		for (int j = 0; j < 3; j++)
//		{
//			Q[j] = point[j] + kSet[i][j];
//		}
//
//		//temp: k在n方向上的投影数值大小
//		double temp = (kSet[i][0] * zNormal[0] + kSet[i][1] * zNormal[1] + kSet[i][2] * zNormal[2]) /
//			sqrt(zNormal[0] * zNormal[0] + zNormal[1] * zNormal[1] + zNormal[2] * zNormal[2]);
//
//		//单位化zNormal
//		double n[3];
//		n[0] = normalized(zNormal)[0];
//		n[1] = normalized(zNormal)[1];
//		n[2] = normalized(zNormal)[2];
//
//		double j[3];
//		for (int k = 0; k < 3; k++)
//		{
//			j[k] = temp * n[k];
//		}
//
//		//R: 投影后的点R
//		double R[3];
//		for (int k = 0; k < 3; k++)
//		{
//			R[k] = Q[k] - j[k];
//		}
//
//		//r: 投影后向量方向
//		double r[3];
//		for (int k = 0; k < 3; k++)
//		{
//			r[k] = R[k] - point[k];
//		}
//
//		//结果写入projNor
//		for (int k = 0; k < 3; k++)
//		{
//			projNor[i][k] = r[k];
//		}
//	}
//
//	return projNor;
//}
//
//int main()
//{
//	double point[3] = { 0.0, 0.0, 0.0 };
//	double zNormal[3] = { 0.0, -0.789316, -0.613987 };
//
//	vector<vector<double>> kSet = { {1.0, 0.0, 0.0},
//									{0.0, 1.0, 0.0} };
//
//	vector<vector<double>> projNor;
//	projNor = projectionNormal(zNormal, kSet, point);
//
//	for (int i = 0; i < 3; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//		std::cout << projNor[i][j] << std::endl;
//		}
//	}
//
//	return 0;
//}



// 读取点云文件

//#include <vtkVertexGlyphFilter.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <sstream>
//
//#include <vtkAutoInit.h> 
//
//int main(int, char* [])
//{
//	std::string filename = "point.txt";
//	std::ifstream filestream(filename.c_str());
//
//	std::string line;
//	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//	while (std::getline(filestream, line))
//	{
//		double x, y, z;
//		std::stringstream linestream;
//		linestream << line;
//		linestream >> x >> y >> z;
//
//		points->InsertNextPoint(x, y, z);
//	}
//
//	filestream.close();
//
//	vtkSmartPointer<vtkPolyData> polyData =
//		vtkSmartPointer<vtkPolyData>::New();
//
//	polyData->SetPoints(points);
//
//	vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	glyphFilter->SetInputData(polyData);
//	glyphFilter->Update();
//
//	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(glyphFilter->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//	renderer->AddActor(actor);
//	renderer->SetBackground(.3, .6, .3);
//
//	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
//	renWin->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	iren->SetRenderWindow(renWin);
//	vtkSmartPointer< vtkInteractorStyleTrackballCamera> style = vtkSmartPointer< vtkInteractorStyleTrackballCamera>::New();
//	iren->SetInteractorStyle(style);
//
//	renWin->SetSize(600, 600);
//	renWin->Render();
//	iren->Start();
//
//	return 0;
//}



//三维一层方向设置

//#include <cmath>
//#include <iostream>
//#include <vector>
//
//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataNormals.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkCellData.h>
//#include <vtkPointData.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkConeSource.h>
//#include <vtkTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkGlyph3D.h>
//#include <vtkGlyphSource2D.h>
//#include <vtkCellCenters.h>
//#include <vtkIdList.h>
//
//#include <algorithm>
//#include <array>
//#include <string>
//#include <fstream>
//
//
//class Vector3
//{
//public:
//	Vector3();
//	Vector3(double X_, double Y_, double Z_)
//	{
//		m_x = X_;
//		m_y = Y_;
//		m_z = Z_;
//	}
//
//	const double x() { return m_x; }
//	const double y() { return m_y; }
//	const double z() { return m_z; }
//
//	Vector3 cross(Vector3& vec)
//	{
//		double X_ = this->y() * vec.z() - this->z() * vec.y();
//		double Y_ = this->z() * vec.x() - this->x() * vec.z();
//		double Z_ = this->x() * vec.y() - this->y() * vec.x();
//
//		return Vector3(X_, Y_, Z_);
//	}
//
//	void normalized()
//	{
//		double w = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
//		m_x /= w;
//		m_y /= w;
//		m_z /= w;
//	}
//
//private:
//	double m_x;
//	double m_y;
//	double m_z;
//};
//
////给定三个点p1 p2 p3，返回由这三个点确定平面的法向n，已归一化
//Vector3 getNormal(vector<double> p1, vector<double> p2, vector<double> p3)
//{
//	Vector3 v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
//	Vector3 v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
//
//	Vector3 n = v1.cross(v2);
//	n.normalized();
//
//	return n;
//}
//
//vector<double> normalized(vector<double>& vector)
//{
//	double L = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
//	vector[0] /= L;
//	vector[1] /= L;
//	vector[2] /= L;
//	return vector;
//}
//
///*zNormal:该单元法向
//* kSet:用户设置的x或y方向
//* point:该单元中心坐标
//* 返回x或y方向渗透率的投影方向
//*/
//vector<double>  projectionNormal(vector<double> zNormal, const vector<double> kSet, const double* point)
//{
//	//分别计算xNor（kSet[0]）和yNor（kSet[1]）
//	double Q[3];
//	for (int i = 0; i < 3; i++)
//	{
//		Q[i] = point[i] + kSet[i];
//	}
//
//	//temp: k在n方向上的投影数值大小
//	double temp = (kSet[0] * zNormal[0] + kSet[1] * zNormal[1] + kSet[2] * zNormal[2]) /
//		sqrt(zNormal[0] * zNormal[0] + zNormal[1] * zNormal[1] + zNormal[2] * zNormal[2]);
//
//	//单位化zNormal
//	double n[3];
//	n[0] = normalized(zNormal)[0];
//	n[1] = normalized(zNormal)[1];
//	n[2] = normalized(zNormal)[2];
//
//	double j[3];
//	for (int i = 0; i < 3; i++)
//	{
//		j[i] = temp * n[i];
//	}
//
//	//R: 投影后的点R
//	double R[3];
//	for (int i = 0; i < 3; i++)
//	{
//		R[i] = Q[i] - j[i];
//	}
//
//	//r: 投影后向量方向
//	double r[3];
//	for (int i = 0; i < 3; i++)
//	{
//		r[i] = R[i] - point[i];
//	}
//
//	//结果写入projNor
//	vector<double> projNor(3);
//	for (int i = 0; i < 3; i++)
//	{
//		projNor[i] = r[i];
//	}
//
//	return projNor;
//}
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	//renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->SetBackground(colors->GetColor3d("Black").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//reader->SetFileName("Test01.vtk");
//	//reader->SetFileName("step12.vtk");
//	//reader->SetFileName("blade.vtk");
//	//reader->SetFileName("plane.vtk");
//	//reader->SetFileName("S.vtk");
//	//reader->SetFileName("qu.vtk");
//	reader->SetFileName("blade2lay.vtk");
//	//reader->SetFileName("bladeGmsh3_vtk.vtk");
//	reader->Update();
//
//	auto unstructuredGrid = reader->GetOutput();
//
//	// 叉乘得到的法向
//	vector<double> calNor;
//
//	// 单元中心坐标
//	vector<double> centrePoint;
//
//	for (int i = 0; i < unstructuredGrid->GetNumberOfCells(); i++)
//	{
//		vector<double> p1, p2, p3, p4, p5, p6, p7, p8;
//
//		//get单元i包含的所有点的id
//		vtkIdList* pointList = unstructuredGrid->GetCell(i)->GetPointIds();
//		int id0 = pointList->GetId(0);
//		int id1 = pointList->GetId(1);
//		int id2 = pointList->GetId(2);
//		int id3 = pointList->GetId(3);
//		int id4 = pointList->GetId(4);
//		int id5 = pointList->GetId(5);
//		int id6 = pointList->GetId(6);
//		int id7 = pointList->GetId(7);
//
//		//给单元i包含的点id排序，arr[]为排序后id编号数组
//		int arr[8] = { id0, id1, id2, id3, id4, id5, id6, id7 };
//		for (int j = 0; j < 8; j++)
//		{
//			int value = arr[j];
//			int position = j;
//			while (position > 0 && arr[position - 1] > value)
//			{
//				arr[position] = arr[position - 1];
//				position--;
//			}
//			arr[position] = value;
//		}
//
//		double* coors1 = unstructuredGrid->GetPoint(arr[0]);
//		p1.push_back(*coors1);
//		p1.push_back(*(coors1 + 1));
//		p1.push_back(*(coors1 + 2));
//
//		double* coors2 = unstructuredGrid->GetPoint(arr[1]);
//		p2.push_back(*coors2);
//		p2.push_back(*(coors2 + 1));
//		p2.push_back(*(coors2 + 2));
//
//		double* coors3 = unstructuredGrid->GetPoint(arr[2]);
//		p3.push_back(*coors3);
//		p3.push_back(*(coors3 + 1));
//		p3.push_back(*(coors3 + 2));
//
//		double* coors4 = unstructuredGrid->GetPoint(arr[3]);
//		p4.push_back(*coors4);
//		p4.push_back(*(coors4 + 1));
//		p4.push_back(*(coors4 + 2));
//
//		double* coors5 = unstructuredGrid->GetPoint(arr[4]);
//		p5.push_back(*coors5);
//		p5.push_back(*(coors5 + 1));
//		p5.push_back(*(coors5 + 2));
//
//		double* coors6 = unstructuredGrid->GetPoint(arr[5]);
//		p6.push_back(*coors6);
//		p6.push_back(*(coors6 + 1));
//		p6.push_back(*(coors6 + 2));
//
//		double* coors7 = unstructuredGrid->GetPoint(arr[6]);
//		p7.push_back(*coors7);
//		p7.push_back(*(coors7 + 1));
//		p7.push_back(*(coors7 + 2));
//
//		double* coors8 = unstructuredGrid->GetPoint(arr[7]);
//		p8.push_back(*coors8);
//		p8.push_back(*(coors8 + 1));
//		p8.push_back(*(coors8 + 2));
//
//		Vector3 cal_nor = getNormal(p1, p2, p3);
//
//		// 把p1 p2 p3平面的法向按顺序存进calNor
//		calNor.push_back(cal_nor.x());
//		calNor.push_back(cal_nor.y());
//		calNor.push_back(cal_nor.z());
//
//		double x = (p1[0] + p2[0] + p3[0] + p4[0] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//		double y = (p1[1] + p2[1] + p3[1] + p4[1] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//		double z = (p1[2] + p2[2] + p3[2] + p4[2] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//
//
//		// 把i单元的中心点按顺序存进centrePoint
//		centrePoint.push_back(x);
//		centrePoint.push_back(y);
//		centrePoint.push_back(z);
//	}
//
//	ofstream ofs2;
//	ofs2.open("calculate_Normal.txt", ios::out);
//	ofs2 << "For points ---------------------: \n";
//
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		ofs2 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	// unstructuredGrid 转换成 polydata
//	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
//		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//	surfaceFilter->SetInputData(unstructuredGrid);
//	surfaceFilter->Update();
//
//	vtkPolyData* polydata = surfaceFilter->GetOutput();
//
//	/*------------------------------数据可视化------------------------------*/
//	// Visualize
//	vtkNew<vtkDataSetMapper> dataSetMapper;
//	dataSetMapper->SetInputData(unstructuredGrid);
//	dataSetMapper->ScalarVisibilityOff();
//
//	vtkNew<vtkProperty> backProp;
//	backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
//	backProp->SetSpecular(.6);
//	backProp->SetSpecularPower(30);
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(dataSetMapper);
//	actor->SetBackfaceProperty(backProp);
//	actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
//	actor->GetProperty()->SetSpecular(.3);
//	actor->GetProperty()->SetSpecularPower(30);
//	actor->GetProperty()->EdgeVisibilityOn();
//	renderer->AddActor(actor);
//	/*-----------------------------------------------------------------*/
//
//	vtkSmartPointer<vtkPolyDataNormals> pdNormals =
//		vtkSmartPointer<vtkPolyDataNormals>::New();
//	pdNormals->SetInputConnection(surfaceFilter->GetOutputPort());
//	pdNormals->ComputeCellNormalsOn();
//	pdNormals->Update();
//
//	vtkPointData* ptData = pdNormals->GetOutput()->GetPointData();
//
//	ofstream ofs;
//	vector<double*> zNormal;
//	vector<double> kx;
//	vector<double> ky;
//	vector<double> kz;
//
//	//两个方向kx和ky
//	vector<double> kxSet{ 1.0, 0.0, 0.0 };
//	vector<double> kySet{ 0.0, 1.0, 0.0 };
//
//	if (ptData)
//	{
//		vtkDataArray* ptNormals = pdNormals->GetOutput()->GetPointData()->GetNormals();
//		if (ptNormals)
//		{
//			//ofs创建normal.txt，写入文件头
//			ofs.open("normal.txt", ios::out);
//			ofs << "For points in every cell: \n";
//			ofs << ptNormals->GetNumberOfTuples() << endl;
//
//			cout << "For points in every cell: \n";
//			cout << ptNormals->GetNumberOfTuples() << endl;
//			for (int i = 0; i < ptNormals->GetNumberOfTuples(); ++i)
//			{
//				double value[3];	//索引为i的法线方向
//				ptNormals->GetTuple(i, value);
//
//				// 输出点上法线方向
//				ofs << i << " " << value[0] << " " << value[1] << " " << value[2] << endl;
//
//				zNormal.push_back(value);	//装入zNormal
//
//				// 是否换成kx、ky、kz ???  是
//				kz.push_back(*value);
//				kz.push_back(*(value + 1));
//				kz.push_back(*(value + 2));
//
//				//	测试zNormal是否正确：正确
//				std::cout << "ZNORMAL COUT::  " << (double)zNormal[i][0] << " " << (double)zNormal[i][1] << " " << (double)zNormal[i][2] << std::endl;
//				cout << "KZ COUT:  " << kz[(0 + 3 * i)] << " " << kz[(1 + 3 * i)] << " " << kz[(2 + 3 * i)] << endl;
//				cout << endl;
//			}
//
//			//测试循环外Kz数值是否正确：正确
//			for (int k = 0; k < kz.size(); k += 3)
//			{
//				cout << "AFTER FOR KZ COUT:  " << kz[k] << " " << kz[k + 1] << " " << kz[k + 2] << endl;
//			}
//		}
//	}
//
//	ofstream ofs1;
//	ofs1.open("normalAndPoint.txt", ios::out);
//	ofs1 << "For points ---------------------: \n";
//
//	// 投影
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		// double* point = pdNormals->GetOutput()->GetPoint(i);
//		double point[3] = { centrePoint[i], centrePoint[i + 1], centrePoint[i + 2] };
//		vector<double> n{ calNor[i], calNor[i + 1], calNor[i + 2] };
//
//		ofs1 << point[0] << " " << point[1] << " " << point[2] << endl;
//
//		vector<double> xProjNor(3);	//x方向上的投影
//
//		xProjNor = projectionNormal(n, kxSet, point);
//		// 计算出的kxSet放入kx中
//		kx.push_back(xProjNor[0]);
//		kx.push_back(xProjNor[1]);
//		kx.push_back(xProjNor[2]);
//
//		vector<double> yProjNor(3);	// y方向上投影
//
//		yProjNor = projectionNormal(n, kySet, point);
//		// 计算出的kySet放入ky中
//		ky.push_back(yProjNor[0]);
//		ky.push_back(yProjNor[1]);
//		ky.push_back(yProjNor[2]);
//	}
//
//	ofs1 << "For Normals  kx  ---------------------: \n";
//	for (int i = 0; i < kx.size(); i += 3)
//	{
//		ofs1 << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  ky  ---------------------: \n";
//	for (int i = 0; i < ky.size(); i += 3)
//	{
//		ofs1 << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  kz  ---------------------: \n";
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		ofs1 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	for (int i = 0; i < kx.size(); i += 3)
//	{
//		ofs << i / 3 << " " << "kx: " << " " << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	for (int i = 0; i < ky.size(); i += 3)
//	{
//		ofs << i / 3 << " " << "ky: " << " " << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//	vtkNew<vtkGlyphSource2D>arrow;//二维箭头
//	arrow->SetGlyphTypeToArrow();
//	arrow->FilledOff();
//
//	vtkSmartPointer<vtkTransform> transform =
//		vtkSmartPointer<vtkTransform>::New();
//	transform->RotateY(180); // make vertex outside
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformF =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	//transformF->SetInputConnection(cone->GetOutputPort());
//	transformF->SetInputConnection(arrow->GetOutputPort());
//	transformF->SetTransform(transform);
//
//	vtkSmartPointer<vtkGlyph3D> glyph =
//		vtkSmartPointer<vtkGlyph3D>::New();
//	glyph->SetInputConnection(pdNormals->GetOutputPort());
//	//glyph->SetSourceConnection(transformF->GetOutputPort()); // source => transform => graph3D
//	glyph->SetSourceConnection(arrow->GetOutputPort()); // source => transform => graph3D
//	glyph->SetVectorModeToUseNormal();
//	glyph->SetScaleModeToScaleByVector();
//	glyph->SetScaleFactor(16);
//
//	vtkSmartPointer<vtkPolyDataMapper> spikeMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	spikeMapper->SetInputConnection(glyph->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> spikeActor = vtkSmartPointer<vtkActor>::New();
//	spikeActor->SetMapper(spikeMapper);
//	spikeActor->GetProperty()->SetColor(0.0, 0.79, 0.34);
//
//	// Add the actors to the renderer, set the background and size
//	//renderer->AddActor(surfaceActor);
//	renderer->AddActor(spikeActor);
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}



//三维多层方向设置  适用于drag拉伸的网格

//#include <cmath>
//#include <iostream>
//#include <vector>
//
//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataNormals.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkCellData.h>
//#include <vtkPointData.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkConeSource.h>
//#include <vtkTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkGlyph3D.h>
//#include <vtkGlyphSource2D.h>
//#include <vtkCellCenters.h>
//#include <vtkIdList.h>
//
//#include <algorithm>
//#include <array>
//#include <string>
//#include <fstream>
//
//
//class Vector3
//{
//public:
//	Vector3();
//	Vector3(double X_, double Y_, double Z_)
//	{
//		m_x = X_;
//		m_y = Y_;
//		m_z = Z_;
//	}
//
//	const double x() { return m_x; }
//	const double y() { return m_y; }
//	const double z() { return m_z; }
//
//	Vector3 cross(Vector3& vec)
//	{
//		double X_ = this->y() * vec.z() - this->z() * vec.y();
//		double Y_ = this->z() * vec.x() - this->x() * vec.z();
//		double Z_ = this->x() * vec.y() - this->y() * vec.x();
//
//		return Vector3(X_, Y_, Z_);
//	}
//
//	void normalized()
//	{
//		double w = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
//		m_x /= w;
//		m_y /= w;
//		m_z /= w;
//	}
//
//private:
//	double m_x;
//	double m_y;
//	double m_z;
//};
//
////给定三个点p1 p2 p3，返回由这三个点确定平面的法向n，已归一化
//Vector3 getNormal(vector<double> p1, vector<double> p2, vector<double> p3)
//{
//	Vector3 v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
//	Vector3 v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
//
//	Vector3 n = v1.cross(v2);
//	n.normalized();
//
//	return n;
//}
//
//vector<double> normalized(vector<double>& vector)
//{
//	double L = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
//	vector[0] /= L;
//	vector[1] /= L;
//	vector[2] /= L;
//	return vector;
//}
//
///*zNormal:该单元法向
//* kSet:用户设置的x或y方向
//* point:该单元中心坐标
//* 返回x或y方向渗透率的投影方向
//*/
//vector<double>  projectionNormal(vector<double> zNormal, const vector<double> kSet, const double* point)
//{
//	//分别计算xNor（kSet[0]）和yNor（kSet[1]）
//	double Q[3];
//	for (int i = 0; i < 3; i++)
//	{
//		Q[i] = point[i] + kSet[i];
//	}
//
//	//temp: k在n方向上的投影数值大小
//	double temp = (kSet[0] * zNormal[0] + kSet[1] * zNormal[1] + kSet[2] * zNormal[2]) /
//		sqrt(zNormal[0] * zNormal[0] + zNormal[1] * zNormal[1] + zNormal[2] * zNormal[2]);
//
//	//单位化zNormal
//	double n[3];
//	n[0] = normalized(zNormal)[0];
//	n[1] = normalized(zNormal)[1];
//	n[2] = normalized(zNormal)[2];
//
//	double j[3];
//	for (int i = 0; i < 3; i++)
//	{
//		j[i] = temp * n[i];
//	}
//
//	//R: 投影后的点R
//	double R[3];
//	for (int i = 0; i < 3; i++)
//	{
//		R[i] = Q[i] - j[i];
//	}
//
//	//r: 投影后向量方向
//	double r[3];
//	for (int i = 0; i < 3; i++)
//	{
//		r[i] = R[i] - point[i];
//	}
//
//	//结果写入projNor
//	vector<double> projNor(3);
//	for (int i = 0; i < 3; i++)
//	{
//		projNor[i] = r[i];
//	}
//
//	return projNor;
//}
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	//renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->SetBackground(colors->GetColor3d("Black").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//reader->SetFileName("Test01.vtk");
//	//reader->SetFileName("step12.vtk");
//	//reader->SetFileName("blade.vtk");
//	//reader->SetFileName("plane.vtk");
//	//reader->SetFileName("S.vtk");
//	//reader->SetFileName("qu.vtk");
//	reader->SetFileName("blade3.vtk");
//	//reader->SetFileName("bladeGmsh3_vtk.vtk");
//	reader->Update();
//
//	auto unstructuredGrid = reader->GetOutput();
//
//	// 叉乘得到的法向
//	vector<double> calNor;
//
//	// 单元中心坐标
//	vector<double> centrePoint;
//
//	for (int i = 0; i < unstructuredGrid->GetNumberOfCells(); i++)
//	{
//		vector<double> p1, p2, p3, p4, p5, p6, p7, p8;
//
//		//get单元i包含的所有点的id
//		vtkIdList* pointList = unstructuredGrid->GetCell(i)->GetPointIds();
//		int id0 = pointList->GetId(0);
//		int id1 = pointList->GetId(1);
//		int id2 = pointList->GetId(2);
//		int id3 = pointList->GetId(3);
//		int id4 = pointList->GetId(4);
//		int id5 = pointList->GetId(5);
//		int id6 = pointList->GetId(6);
//		int id7 = pointList->GetId(7);
//
//		//给单元i包含的点id排序，arr[]为排序后id编号数组
//		int arr[8] = { id0, id1, id2, id3, id4, id5, id6, id7 };
//		for (int j = 0; j < 8; j++)
//		{
//			int value = arr[j];
//			int position = j;
//			while (position > 0 && arr[position - 1] > value)
//			{
//				arr[position] = arr[position - 1];
//				position--;
//			}
//			arr[position] = value;
//		}
//
//		double* coors1 = unstructuredGrid->GetPoint(arr[0]);
//		p1.push_back(*coors1);
//		p1.push_back(*(coors1 + 1));
//		p1.push_back(*(coors1 + 2));
//
//		double* coors2 = unstructuredGrid->GetPoint(arr[1]);
//		p2.push_back(*coors2);
//		p2.push_back(*(coors2 + 1));
//		p2.push_back(*(coors2 + 2));
//
//		double* coors3 = unstructuredGrid->GetPoint(arr[2]);
//		p3.push_back(*coors3);
//		p3.push_back(*(coors3 + 1));
//		p3.push_back(*(coors3 + 2));
//
//		double* coors4 = unstructuredGrid->GetPoint(arr[3]);
//		p4.push_back(*coors4);
//		p4.push_back(*(coors4 + 1));
//		p4.push_back(*(coors4 + 2));
//
//		double* coors5 = unstructuredGrid->GetPoint(arr[4]);
//		p5.push_back(*coors5);
//		p5.push_back(*(coors5 + 1));
//		p5.push_back(*(coors5 + 2));
//
//		double* coors6 = unstructuredGrid->GetPoint(arr[5]);
//		p6.push_back(*coors6);
//		p6.push_back(*(coors6 + 1));
//		p6.push_back(*(coors6 + 2));
//
//		double* coors7 = unstructuredGrid->GetPoint(arr[6]);
//		p7.push_back(*coors7);
//		p7.push_back(*(coors7 + 1));
//		p7.push_back(*(coors7 + 2));
//
//		double* coors8 = unstructuredGrid->GetPoint(arr[7]);
//		p8.push_back(*coors8);
//		p8.push_back(*(coors8 + 1));
//		p8.push_back(*(coors8 + 2));
//
//		Vector3 cal_nor = getNormal(p1, p2, p3);
//
//		// 把p1 p2 p3平面的法向按顺序存进calNor
//		calNor.push_back(cal_nor.x());
//		calNor.push_back(cal_nor.y());
//		calNor.push_back(cal_nor.z());
//
//		double x = (p1[0] + p2[0] + p3[0] + p4[0] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//		double y = (p1[1] + p2[1] + p3[1] + p4[1] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//		double z = (p1[2] + p2[2] + p3[2] + p4[2] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//
//
//		// 把i单元的中心点按顺序存进centrePoint
//		centrePoint.push_back(x);
//		centrePoint.push_back(y);
//		centrePoint.push_back(z);
//	}
//
//	ofstream ofs2;
//	ofs2.open("calculate_Normal.txt", ios::out);
//	ofs2 << "For points ---------------------: \n";
//
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		ofs2 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	// unstructuredGrid 转换成 polydata
//	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
//		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//	surfaceFilter->SetInputData(unstructuredGrid);
//	surfaceFilter->Update();
//
//	vtkPolyData* polydata = surfaceFilter->GetOutput();
//
//	/*------------------------------数据可视化------------------------------*/
//	// Visualize
//	vtkNew<vtkDataSetMapper> dataSetMapper;
//	dataSetMapper->SetInputData(unstructuredGrid);
//	dataSetMapper->ScalarVisibilityOff();
//
//	vtkNew<vtkProperty> backProp;
//	backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
//	backProp->SetSpecular(.6);
//	backProp->SetSpecularPower(30);
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(dataSetMapper);
//	actor->SetBackfaceProperty(backProp);
//	actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
//	actor->GetProperty()->SetSpecular(.3);
//	actor->GetProperty()->SetSpecularPower(30);
//	actor->GetProperty()->EdgeVisibilityOn();
//	renderer->AddActor(actor);
//	/*-----------------------------------------------------------------*/
//
//	vtkSmartPointer<vtkPolyDataNormals> pdNormals =
//		vtkSmartPointer<vtkPolyDataNormals>::New();
//	pdNormals->SetInputConnection(surfaceFilter->GetOutputPort());
//	pdNormals->ComputeCellNormalsOn();
//	pdNormals->Update();
//
//	vtkPointData* ptData = pdNormals->GetOutput()->GetPointData();
//
//	ofstream ofs;
//	vector<double*> zNormal;
//	vector<double> kx;
//	vector<double> ky;
//	vector<double> kz;
//
//	//两个方向kx和ky
//	vector<double> kxSet{ 1.0, 0.0, 0.0 };
//	vector<double> kySet{ 0.0, 1.0, 0.0 };
//
//	if (ptData)
//	{
//		vtkDataArray* ptNormals = pdNormals->GetOutput()->GetPointData()->GetNormals();
//		if (ptNormals)
//		{
//			//ofs创建normal.txt，写入文件头
//			ofs.open("normal.txt", ios::out);
//			ofs << "For points in every cell: \n";
//			ofs << ptNormals->GetNumberOfTuples() << endl;
//
//			cout << "For points in every cell: \n";
//			cout << ptNormals->GetNumberOfTuples() << endl;
//			for (int i = 0; i < ptNormals->GetNumberOfTuples(); ++i)
//			{
//				double value[3];	//索引为i的法线方向
//				ptNormals->GetTuple(i, value);
//
//				// 输出点上法线方向
//				ofs << i << " " << value[0] << " " << value[1] << " " << value[2] << endl;
//
//				zNormal.push_back(value);	//装入zNormal
//
//				// 是否换成kx、ky、kz ???  是
//				kz.push_back(*value);
//				kz.push_back(*(value + 1));
//				kz.push_back(*(value + 2));
//
//				//	测试zNormal是否正确：正确
//				std::cout << "ZNORMAL COUT::  " << (double)zNormal[i][0] << " " << (double)zNormal[i][1] << " " << (double)zNormal[i][2] << std::endl;
//				cout << "KZ COUT:  " << kz[(0 + 3 * i)] << " " << kz[(1 + 3 * i)] << " " << kz[(2 + 3 * i)] << endl;
//				cout << endl;
//			}
//
//			//测试循环外Kz数值是否正确：正确
//			for (int k = 0; k < kz.size(); k += 3)
//			{
//				cout << "AFTER FOR KZ COUT:  " << kz[k] << " " << kz[k + 1] << " " << kz[k + 2] << endl;
//			}
//		}
//	}
//
//	ofstream ofs1;
//	ofs1.open("normalAndPoint.txt", ios::out);
//	ofs1 << "For points ---------------------: \n";
//
//	// 投影
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		// double* point = pdNormals->GetOutput()->GetPoint(i);
//		double point[3] = { centrePoint[i], centrePoint[i + 1], centrePoint[i + 2] };
//		vector<double> n{ calNor[i], calNor[i + 1], calNor[i + 2] };
//
//		ofs1 << point[0] << " " << point[1] << " " << point[2] << endl;
//
//		vector<double> xProjNor(3);	//x方向上的投影
//
//		xProjNor = projectionNormal(n, kxSet, point);
//		// 计算出的kxSet放入kx中
//		kx.push_back(xProjNor[0]);
//		kx.push_back(xProjNor[1]);
//		kx.push_back(xProjNor[2]);
//
//		vector<double> yProjNor(3);	// y方向上投影
//
//		yProjNor = projectionNormal(n, kySet, point);
//		// 计算出的kySet放入ky中
//		ky.push_back(yProjNor[0]);
//		ky.push_back(yProjNor[1]);
//		ky.push_back(yProjNor[2]);
//	}
//
//	ofs1 << "For Normals  kx  ---------------------: \n";
//	for (int i = 0; i < kx.size(); i += 3)
//	{
//		ofs1 << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  ky  ---------------------: \n";
//	for (int i = 0; i < ky.size(); i += 3)
//	{
//		ofs1 << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  kz  ---------------------: \n";
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		ofs1 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	for (int i = 0; i < kx.size(); i += 3)
//	{
//		ofs << i / 3 << " " << "kx: " << " " << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	for (int i = 0; i < ky.size(); i += 3)
//	{
//		ofs << i / 3 << " " << "ky: " << " " << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//	vtkNew<vtkGlyphSource2D>arrow;//二维箭头
//	arrow->SetGlyphTypeToArrow();
//	arrow->FilledOff();
//
//	vtkSmartPointer<vtkTransform> transform =
//		vtkSmartPointer<vtkTransform>::New();
//	transform->RotateY(180); // make vertex outside
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformF =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	//transformF->SetInputConnection(cone->GetOutputPort());
//	transformF->SetInputConnection(arrow->GetOutputPort());
//	transformF->SetTransform(transform);
//
//	vtkSmartPointer<vtkGlyph3D> glyph =
//		vtkSmartPointer<vtkGlyph3D>::New();
//	glyph->SetInputConnection(pdNormals->GetOutputPort());
//	//glyph->SetSourceConnection(transformF->GetOutputPort()); // source => transform => graph3D
//	glyph->SetSourceConnection(arrow->GetOutputPort()); // source => transform => graph3D
//	glyph->SetVectorModeToUseNormal();
//	glyph->SetScaleModeToScaleByVector();
//	glyph->SetScaleFactor(16);
//
//	vtkSmartPointer<vtkPolyDataMapper> spikeMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	spikeMapper->SetInputConnection(glyph->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> spikeActor = vtkSmartPointer<vtkActor>::New();
//	spikeActor->SetMapper(spikeMapper);
//	spikeActor->GetProperty()->SetColor(0.0, 0.79, 0.34);
//
//	// Add the actors to the renderer, set the background and size
//	//renderer->AddActor(surfaceActor);
//	renderer->AddActor(spikeActor);
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}



//三维多层方向设置  偏置网格

//#include <cmath>
//#include <iostream>
//#include <vector>
//
//#include <vtkAppendFilter.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkDataSetMapper.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataNormals.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkCellData.h>
//#include <vtkPointData.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkConeSource.h>
//#include <vtkTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkGlyph3D.h>
//#include <vtkGlyphSource2D.h>
//#include <vtkCellCenters.h>
//#include <vtkIdList.h>
//
//#include <algorithm>
//#include <array>
//#include <string>
//#include <fstream>
//
//
//class Vector3
//{
//public:
//	Vector3();
//	Vector3(double X_, double Y_, double Z_)
//	{
//		m_x = X_;
//		m_y = Y_;
//		m_z = Z_;
//	}
//
//	const double x() { return m_x; }
//	const double y() { return m_y; }
//	const double z() { return m_z; }
//
//	Vector3 cross(Vector3& vec)
//	{
//		double X_ = this->y() * vec.z() - this->z() * vec.y();
//		double Y_ = this->z() * vec.x() - this->x() * vec.z();
//		double Z_ = this->x() * vec.y() - this->y() * vec.x();
//
//		return Vector3(X_, Y_, Z_);
//	}
//
//	void normalized()
//	{
//		double w = sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
//		m_x /= w;
//		m_y /= w;
//		m_z /= w;
//	}
//
//private:
//	double m_x;
//	double m_y;
//	double m_z;
//};
//
////给定三个点p1 p2 p3，返回由这三个点确定平面的法向n，已归一化
//Vector3 getNormal(vector<double> p1, vector<double> p2, vector<double> p3)
//{
//	Vector3 v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
//	Vector3 v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
//
//	Vector3 n = v1.cross(v2);
//	n.normalized();
//
//	return n;
//}
//
//vector<double> normalized(vector<double>& vector)
//{
//	double L = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
//	vector[0] /= L;
//	vector[1] /= L;
//	vector[2] /= L;
//	return vector;
//}
//
///*zNormal:该单元法向
//* kSet:用户设置的x或y方向
//* point:该单元中心坐标
//* 返回x或y方向渗透率的投影方向
//*/
//vector<double>  projectionNormal(vector<double> zNormal, const vector<double> kSet, const double* point)
//{
//	//分别计算xNor（kSet[0]）和yNor（kSet[1]）
//	double Q[3];
//	for (int i = 0; i < 3; i++)
//	{
//		Q[i] = point[i] + kSet[i];
//	}
//
//	//temp: k在n方向上的投影数值大小
//	double temp = (kSet[0] * zNormal[0] + kSet[1] * zNormal[1] + kSet[2] * zNormal[2]) /
//		sqrt(zNormal[0] * zNormal[0] + zNormal[1] * zNormal[1] + zNormal[2] * zNormal[2]);
//
//	//单位化zNormal
//	double n[3];
//	n[0] = normalized(zNormal)[0];
//	n[1] = normalized(zNormal)[1];
//	n[2] = normalized(zNormal)[2];
//
//	double j[3];
//	for (int i = 0; i < 3; i++)
//	{
//		j[i] = temp * n[i];
//	}
//
//	//R: 投影后的点R
//	double R[3];
//	for (int i = 0; i < 3; i++)
//	{
//		R[i] = Q[i] - j[i];
//	}
//
//	//r: 投影后向量方向
//	double r[3];
//	for (int i = 0; i < 3; i++)
//	{
//		r[i] = R[i] - point[i];
//	}
//
//	//结果写入projNor
//	vector<double> projNor(3);
//	for (int i = 0; i < 3; i++)
//	{
//		projNor[i] = r[i];
//	}
//
//	return projNor;
//}
//
//namespace {
//	vtkSmartPointer<vtkUnstructuredGrid>
//		ReadUnstructuredGrid(std::string const& fileName);
//}
//
//static int N = 4; // N为厚度方向上网格层数
//
//int main(int argc, char* argv[])
//{
//	// Vis Pipeline
//	vtkNew<vtkNamedColors> colors;
//
//	vtkNew<vtkRenderer> renderer;
//
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->SetSize(640, 480);
//	renderWindow->AddRenderer(renderer);
//
//	vtkNew<vtkRenderWindowInteractor> interactor;
//	interactor->SetRenderWindow(renderWindow);
//
//	//renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
//	renderer->SetBackground(colors->GetColor3d("Black").GetData());
//	renderer->UseHiddenLineRemovalOn();
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	//reader->SetFileName("Test01.vtk");
//	//reader->SetFileName("step12.vtk");
//	//reader->SetFileName("blade.vtk");
//	//reader->SetFileName("plane.vtk");
//	//reader->SetFileName("S.vtk");
//	//reader->SetFileName("qu.vtk");
//	reader->SetFileName("wave1.vtk");
//	//reader->SetFileName("L2.vtk");
//	//reader->SetFileName("bladeGmsh3_vtk.vtk");
//	reader->Update();
//
//	auto unstructuredGrid = reader->GetOutput();
//
//	// 叉乘得到的法向
//	vector<double> calNor;
//
//	// 单元中心坐标
//	vector<double> centrePoint;
//
//	// 按文件内单元所包含节点编号，用每个单元的前三个节点来计算法向
//	// drag输出后的文件单元包含节点的编号顺序为按拉伸顺序从上往下排序，所以前四个节点id一定比后四个节点id小
//	// linear drag按照单元id依次编号，顺序为第一层2D单元，然后在某个2D单元往下拉伸N个并按序编号，再下一个2D单元，重复， \
//	// \ 所以每个单元的后四个节点id为按偏置方向往下的下一个网格的前四个节点id
//	for (int i = 0; i < unstructuredGrid->GetNumberOfCells(); i++)
//	{
//		vector<double> p1, p2, p3, p4, p5, p6, p7, p8;
//
//		//get单元i包含的所有点的id
//		vtkIdList* pointList = unstructuredGrid->GetCell(i)->GetPointIds();
//		int id0 = pointList->GetId(0);
//		int id1 = pointList->GetId(1);
//		int id2 = pointList->GetId(2);
//		int id3 = pointList->GetId(3);
//		int id4 = pointList->GetId(4);
//		int id5 = pointList->GetId(5);
//		int id6 = pointList->GetId(6);
//		int id7 = pointList->GetId(7);
//
//		double* coors1 = unstructuredGrid->GetPoint(id0);
//		p1.push_back(*coors1);
//		p1.push_back(*(coors1 + 1));
//		p1.push_back(*(coors1 + 2));
//
//		double* coors2 = unstructuredGrid->GetPoint(id1);
//		p2.push_back(*coors2);
//		p2.push_back(*(coors2 + 1));
//		p2.push_back(*(coors2 + 2));
//
//		double* coors3 = unstructuredGrid->GetPoint(id2);
//		p3.push_back(*coors3);
//		p3.push_back(*(coors3 + 1));
//		p3.push_back(*(coors3 + 2));
//
//		double* coors4 = unstructuredGrid->GetPoint(id3);
//		p4.push_back(*coors4);
//		p4.push_back(*(coors4 + 1));
//		p4.push_back(*(coors4 + 2));
//
//		double* coors5 = unstructuredGrid->GetPoint(id4);
//		p5.push_back(*coors5);
//		p5.push_back(*(coors5 + 1));
//		p5.push_back(*(coors5 + 2));
//
//		double* coors6 = unstructuredGrid->GetPoint(id5);
//		p6.push_back(*coors6);
//		p6.push_back(*(coors6 + 1));
//		p6.push_back(*(coors6 + 2));
//
//		double* coors7 = unstructuredGrid->GetPoint(id6);
//		p7.push_back(*coors7);
//		p7.push_back(*(coors7 + 1));
//		p7.push_back(*(coors7 + 2));
//
//		double* coors8 = unstructuredGrid->GetPoint(id7);
//		p8.push_back(*coors8);
//		p8.push_back(*(coors8 + 1));
//		p8.push_back(*(coors8 + 2));
//
//		Vector3 cal_nor = getNormal(p1, p2, p3);
//
//		// 把p1 p2 p3平面的法向按顺序存进calNor
//		calNor.push_back(cal_nor.x());
//		calNor.push_back(cal_nor.y());
//		calNor.push_back(cal_nor.z());
//
//		double x = (p1[0] + p2[0] + p3[0] + p4[0] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//		double y = (p1[1] + p2[1] + p3[1] + p4[1] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//		double z = (p1[2] + p2[2] + p3[2] + p4[2] + p5[0] + p6[0] + p7[0] + p8[0]) / 8;
//
//		// 把i单元的中心点按顺序存进centrePoint
//		centrePoint.push_back(x);
//		centrePoint.push_back(y);
//		centrePoint.push_back(z);
//	}
//
//	// ofs2保存第一层单元法向
//	ofstream ofs2;
//	ofs2.open("calculate_Normal_firstLay.txt", ios::out);
//	ofs2 << "For first lay ---------------------: \n";
//
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		ofs2 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	/*------------------------------数据可视化------------------------------*/
//	vtkNew<vtkDataSetMapper> dataSetMapper;
//	dataSetMapper->SetInputData(unstructuredGrid);
//	dataSetMapper->ScalarVisibilityOff();
//
//	vtkNew<vtkProperty> backProp;
//	backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
//	backProp->SetSpecular(.6);
//	backProp->SetSpecularPower(30);
//
//	vtkNew<vtkActor> actor;
//	actor->SetMapper(dataSetMapper);
//	actor->SetBackfaceProperty(backProp);
//	actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
//	actor->GetProperty()->SetSpecular(.3);
//	actor->GetProperty()->SetSpecularPower(30);
//	actor->GetProperty()->EdgeVisibilityOn();
//	renderer->AddActor(actor);
//
//	// 表面法线可视化数据集
//	// unstructuredGrid 转换成 polydata
//	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
//		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//	surfaceFilter->SetInputData(unstructuredGrid);
//	surfaceFilter->Update();
//	vtkPolyData* polydata = surfaceFilter->GetOutput();
//
//	vtkSmartPointer<vtkPolyDataNormals> pdNormals =
//		vtkSmartPointer<vtkPolyDataNormals>::New();
//	pdNormals->SetInputConnection(surfaceFilter->GetOutputPort());
//	pdNormals->ComputeCellNormalsOn();
//	pdNormals->Update();
//
//	vtkPointData* ptData = pdNormals->GetOutput()->GetPointData();
//	/*-----------------------------------------------------------------*/
//
//	vector<double> kx;
//	vector<double> ky;
//
//	//两个方向渗透率kx和ky
//	vector<double> kxSet{ 1.0, 0.0, 0.0 };
//	vector<double> kySet{ 0.0, 1.0, 0.0 };
//
//	//渗透率方向和基准轴夹角
//	double theta = 30;
//	theta = theta / 180 * 3.1415926;
//
//
//	ofstream ofs1;
//	ofs1.open("normalAndPoint.txt", ios::out);
//	ofs1 << "For points ---------------------: \n";
//
//	// 投影
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		double point[3] = { centrePoint[i], centrePoint[i + 1], centrePoint[i + 2] };
//		vector<double> n{ calNor[i], calNor[i + 1], calNor[i + 2] };
//
//		vector<double> xProjNor(3);	//x方向上的投影
//
//		xProjNor = projectionNormal(n, kxSet, point);
//		// 计算出的kxSet放入kx中
//		//kx.push_back(xProjNor[0]);
//		//kx.push_back(xProjNor[1]);
//		//kx.push_back(xProjNor[2]);
//
//		vector<double> yProjNor(3);	// y方向上投影
//
//		yProjNor = projectionNormal(n, kySet, point);
//		// 计算出的kySet放入ky中
//		//ky.push_back(yProjNor[0]);
//		//ky.push_back(yProjNor[1]);
//		//ky.push_back(yProjNor[2]);
//
//		double R[3][3];
//		R[0][0] = cos(theta) + n[0] * n[0] * (1 - cos(theta));
//		R[0][1] = n[0] * n[1] * (1 - cos(theta)) - n[2] * sin(theta);
//		R[0][2] = n[0] * n[2] * (1 - cos(theta)) + n[1] * sin(theta);
//		R[1][0] = n[1] * n[0] * (1 - cos(theta)) + n[2] * sin(theta);
//		R[1][1] = cos(theta) + n[1] * n[1] * (1 - cos(theta));
//		R[1][2] = n[1] * n[2] * (1 - cos(theta)) - n[0] * sin(theta);
//		R[2][0] = n[2] * n[0] * (1 - cos(theta)) - n[1] * sin(theta);
//		R[2][1] = n[2] * n[1] * (1 - cos(theta)) + n[0] * sin(theta);
//		R[2][2] = cos(theta) + n[2] * n[2] * (1 - cos(theta));
//
//		vector<double> x(3);
//		vector<double> y(3);
//
//		x[0] = R[0][0] * xProjNor[0] + R[0][1] * xProjNor[1] + R[0][2] * xProjNor[2];
//		x[1] = R[1][0] * xProjNor[0] + R[1][1] * xProjNor[1] + R[1][2] * xProjNor[2];
//		x[2] = R[2][0] * xProjNor[0] + R[2][1] * xProjNor[1] + R[2][2] * xProjNor[2];
//
//		y[0] = R[0][0] * yProjNor[0] + R[0][1] * yProjNor[1] + R[0][2] * yProjNor[2];
//		y[1] = R[1][0] * yProjNor[0] + R[1][1] * yProjNor[1] + R[1][2] * yProjNor[2];
//		y[2] = R[2][0] * yProjNor[0] + R[2][1] * yProjNor[1] + R[2][2] * yProjNor[2];
//
//		kx.push_back(x[0]);
//		kx.push_back(x[1]);
//		kx.push_back(x[2]);
//
//		ky.push_back(y[0]);
//		ky.push_back(y[1]);
//		ky.push_back(y[2]);
//	}
//
//	ofs1 << "For Normals  kx  ---------------------: \n";
//	for (int i = 0; i < kx.size(); i += 3)
//	{
//		ofs1 << kx[i] << " " << kx[i + 1] << " " << kx[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  ky  ---------------------: \n";
//	for (int i = 0; i < ky.size(); i += 3)
//	{
//		ofs1 << ky[i] << " " << ky[i + 1] << " " << ky[i + 2] << endl;
//	}
//
//	ofs1 << "For Normals  kz  ---------------------: \n";
//	for (int i = 0; i < calNor.size(); i += 3)
//	{
//		ofs1 << calNor[i] << " " << calNor[i + 1] << " " << calNor[i + 2] << endl;
//	}
//
//	vtkNew<vtkGlyphSource2D>arrow;//二维箭头
//	arrow->SetGlyphTypeToArrow();
//	arrow->FilledOff();
//
//	vtkSmartPointer<vtkTransform> transform =
//		vtkSmartPointer<vtkTransform>::New();
//	transform->RotateY(180); // make vertex outside
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformF =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	//transformF->SetInputConnection(cone->GetOutputPort());
//	transformF->SetInputConnection(arrow->GetOutputPort());
//	transformF->SetTransform(transform);
//
//	vtkSmartPointer<vtkGlyph3D> glyph =
//		vtkSmartPointer<vtkGlyph3D>::New();
//	glyph->SetInputConnection(pdNormals->GetOutputPort());
//	//glyph->SetSourceConnection(transformF->GetOutputPort()); // source => transform => graph3D
//	glyph->SetSourceConnection(arrow->GetOutputPort()); // source => transform => graph3D
//	glyph->SetVectorModeToUseNormal();
//	glyph->SetScaleModeToScaleByVector();
//	glyph->SetScaleFactor(16);
//
//	vtkSmartPointer<vtkPolyDataMapper> spikeMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	spikeMapper->SetInputConnection(glyph->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> spikeActor = vtkSmartPointer<vtkActor>::New();
//	spikeActor->SetMapper(spikeMapper);
//	spikeActor->GetProperty()->SetColor(0.0, 0.79, 0.34);
//
//	// Add the actors to the renderer, set the background and size
//	//renderer->AddActor(surfaceActor);
//	renderer->AddActor(spikeActor);
//	renderer->GetActiveCamera()->Azimuth(45);
//	renderer->GetActiveCamera()->Elevation(45);
//	renderer->ResetCamera();
//	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
//	renderWindow->Render();
//	interactor->Start();
//
//	return EXIT_SUCCESS;
//}




//框选

#include <vtkAppendFilter.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <algorithm>
#include <array>
#include <string>

#include "vtkPolyDataReader.h"
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkProperty.h>
#include <vtkObjectFactory.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include "vtkDecimatePro.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkAppendFilter.h"
#include "vtkDataSetCollection.h"
#include "vtkCollection.h"
#include "vtkType.h"
#include "vtkDataSet.h"
#include "vtkAreaPicker.h"
#include "vtkInteractorStyleRubberBandPick.h"
#include "vtkPlanes.h"
#include "vtkExtractGeometry.h"
#include <vtkPicker.h>
#include <vtkExtractSelectedFrustum.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkPlaneSource.h>

namespace {
	vtkSmartPointer<vtkUnstructuredGrid>
		ReadUnstructuredGrid(std::string const& fileName);
}

class InteractorStyle : public vtkInteractorStyleRubberBandPick
{
public:
	static InteractorStyle* New();
	vtkTypeMacro(InteractorStyle, vtkInteractorStyleRubberBandPick);

	InteractorStyle()
	{
		selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
		selectedActor = vtkSmartPointer<vtkActor>::New();
	}

	virtual void OnLeftButtonUp()
	{
		// Forward events
		vtkNew<vtkNamedColors> colors;
		vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

		vtkPlanes* frustum = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetFrustum();
		vtkPoints* points = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetClipPoints();
		cout << points->GetPoint(0)[0] << " " << points->GetPoint(0)[1] << " " << points->GetPoint(0)[2] << endl;
		cout << points->GetPoint(1)[0] << " " << points->GetPoint(1)[1] << " " << points->GetPoint(1)[2] << endl;
		cout << points->GetPoint(2)[0] << " " << points->GetPoint(2)[1] << " " << points->GetPoint(2)[2] << endl;
		cout << points->GetPoint(3)[0] << " " << points->GetPoint(3)[1] << " " << points->GetPoint(3)[2] << endl;
		cout << points->GetPoint(4)[0] << " " << points->GetPoint(4)[1] << " " << points->GetPoint(4)[2] << endl;
		cout << points->GetPoint(5)[0] << " " << points->GetPoint(5)[1] << " " << points->GetPoint(5)[2] << endl;
		cout << points->GetPoint(6)[0] << " " << points->GetPoint(6)[1] << " " << points->GetPoint(6)[2] << endl;
		cout << points->GetPoint(7)[0] << " " << points->GetPoint(7)[1] << " " << points->GetPoint(7)[2] << endl;
		cout << endl;

		vtkSmartPointer<vtkExtractGeometry> extractGeometry =
			vtkSmartPointer<vtkExtractGeometry>::New();
		extractGeometry->SetImplicitFunction(frustum);

		extractGeometry->SetInputData(this->Data);

		extractGeometry->Update();
		this->selectedMapper->SetInputConnection(extractGeometry->GetOutputPort());
		this->selectedMapper->ScalarVisibilityOff();

		this->selectedActor->SetMapper(selectedMapper);

		this->selectedActor->GetProperty()->EdgeVisibilityOn();
		this->selectedActor->GetProperty()->SetColor(
			colors->GetColor3d("Tomato").GetData());

		// this->selectedActor->GetProperty()->SetLineWidth(3);

		//this->CurrentRenderer->AddActor(SelectedActor);
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
		this->GetInteractor()->GetRenderWindow()->Render();
	}

	vtkSmartPointer<vtkUnstructuredGrid> Data;
	vtkSmartPointer<vtkDataSetMapper> selectedMapper;
	vtkSmartPointer<vtkActor> selectedActor;

};


vtkStandardNewMacro(InteractorStyle);


int main(int argc, char* argv[])
{
	// Vis Pipeline
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(640, 480);
	renderWindow->AddRenderer(renderer);

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);

	renderer->SetBackground(colors->GetColor3d("PaleTurquoise").GetData());
	renderer->UseHiddenLineRemovalOn();

	vtkNew<vtkUnstructuredGridReader> reader;
	reader->SetFileName("wave1.vtk");
	reader->Update();

	auto unstructuredGrid = reader->GetOutput();

	// Visualize
	vtkNew<vtkDataSetMapper> mapper;
	mapper->SetInputData(unstructuredGrid);
	mapper->ScalarVisibilityOff();

	vtkNew<vtkProperty> backProp;
	backProp->SetDiffuseColor(colors->GetColor3d("SeaGreen").GetData());
	backProp->SetSpecular(.6);
	backProp->SetSpecularPower(30);

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->SetBackfaceProperty(backProp);
	actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("SeaGreen").GetData());
	actor->GetProperty()->SetSpecular(.3);
	actor->GetProperty()->SetSpecularPower(30);
	actor->GetProperty()->EdgeVisibilityOn();
	renderer->AddActor(actor);
	renderer->GetActiveCamera()->Azimuth(45);
	renderer->GetActiveCamera()->Elevation(45);
	renderer->ResetCamera();

	vtkSmartPointer<vtkAreaPicker> areaPicker = vtkSmartPointer<vtkAreaPicker>::New();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->SetPicker(areaPicker);
	renderWindowInteractor->Initialize();

	// Set the custom stype to use for interaction.
	vtkSmartPointer<InteractorStyle> style =
		vtkSmartPointer<InteractorStyle>::New();
	//style->SetDefaultRenderer(renderer);
	style->Data = unstructuredGrid;

	renderWindowInteractor->SetInteractorStyle(style);

	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
	renderWindow->Render();
	interactor->Start();

	return EXIT_SUCCESS;
}




//cellPicker

//#include <vtkActor.h>
//#include <vtkCellArray.h>
//#include <vtkCellPicker.h>
//#include <vtkCommand.h>
//#include <vtkDataSetMapper.h>
//#include <vtkExtractSelection.h>
//#include <vtkIdTypeArray.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkObjectFactory.h>
//#include <vtkPlaneSource.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkRendererCollection.h>
//#include <vtkSelection.h>
//#include <vtkSelectionNode.h>
//#include <vtkSmartPointer.h>
//#include <vtkTriangleFilter.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredDataReader.h>
//
//namespace {
//
//	// Catch mouse events
//	class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
//	{
//	public:
//		static MouseInteractorStyle* New();
//
//		MouseInteractorStyle()
//		{
//			selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
//			selectedActor = vtkSmartPointer<vtkActor>::New();
//		}
//
//		virtual void OnLeftButtonDown() override
//		{
//			vtkNew<vtkNamedColors> colors;
//
//			// Get the location of the click (in window coordinates)
//			int* pos = this->GetInteractor()->GetEventPosition();
//
//			vtkNew<vtkCellPicker> picker;
//			picker->SetTolerance(0.0005);
//
//			// Pick from this location.
//			picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
//
//			double* worldPosition = picker->GetPickPosition();
//			std::cout << "Cell id is: " << picker->GetCellId() << std::endl;
//			std::cout << "Pick Screen is: (" << pos[0] << ", "
//				<< pos[1] << ", " << 0 << ")" << endl;
//
//			if (picker->GetCellId() != -1)
//			{
//
//				std::cout << "Pick position is: (" << worldPosition[0] << ", "
//					<< worldPosition[1] << ", " << worldPosition[2] << ")" << endl;
//
//				vtkNew<vtkIdTypeArray> ids;
//				ids->SetNumberOfComponents(1);
//				ids->InsertNextValue(picker->GetCellId());
//
//				vtkNew<vtkSelectionNode> selectionNode;
//				selectionNode->SetFieldType(vtkSelectionNode::CELL);
//				selectionNode->SetContentType(vtkSelectionNode::INDICES);
//				selectionNode->SetSelectionList(ids);
//
//				vtkNew<vtkSelection> selection;
//				selection->AddNode(selectionNode);
//
//				vtkNew<vtkExtractSelection> extractSelection;
//				extractSelection->SetInputData(0, this->Data);
//				extractSelection->SetInputData(1, selection);
//				extractSelection->Update();
//
//				// In selection
//				vtkNew<vtkUnstructuredGrid> selected;
//				selected->ShallowCopy(extractSelection->GetOutput());
//
//				std::cout << "Number of points in the selection: "
//					<< selected->GetNumberOfPoints() << std::endl;
//				std::cout << "Number of cells in the selection : "
//					<< selected->GetNumberOfCells() << std::endl;
//				selectedMapper->SetInputData(selected);
//				selectedActor->SetMapper(selectedMapper);
//				selectedActor->GetProperty()->EdgeVisibilityOn();
//				selectedActor->GetProperty()->SetColor(
//					colors->GetColor3d("Tomato").GetData());
//
//				selectedActor->GetProperty()->SetLineWidth(3);
//
//				this->Interactor->GetRenderWindow()
//					->GetRenderers()
//					->GetFirstRenderer()
//					->AddActor(selectedActor);
//			}
//			// Forward events
//			vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
//		}
//
//		//vtkSmartPointer<vtkPolyData> Data;
//		vtkSmartPointer<vtkUnstructuredGrid> Data;
//		vtkSmartPointer<vtkDataSetMapper> selectedMapper;
//		vtkSmartPointer<vtkActor> selectedActor;
//	};
//
//	vtkStandardNewMacro(MouseInteractorStyle);
//
//} // namespace
//
//int main(int, char* [])
//{
//	vtkNew<vtkNamedColors> colors;
//
//	//vtkNew<vtkPlaneSource> planeSource;
//	//planeSource->Update();
//
//	//vtkNew<vtkTriangleFilter> triangleFilter;
//	//triangleFilter->SetInputConnection(planeSource->GetOutputPort());
//	//triangleFilter->Update();
//
//	//vtkNew<vtkPolyDataMapper> mapper;
//	//mapper->SetInputConnection(triangleFilter->GetOutputPort());
//
//	vtkNew<vtkUnstructuredGridReader> reader;
//	reader->SetFileName("wave1.vtk");
//	reader->Update();
//
//	auto unstructuredGrid = reader->GetOutput();
//
//	// Visualize
//	vtkNew<vtkDataSetMapper> mapper;
//	mapper->SetInputData(unstructuredGrid);
//	mapper->ScalarVisibilityOff();
//
//
//	vtkNew<vtkActor> actor;
//	actor->GetProperty()->SetColor(colors->GetColor3d("SeaGreen").GetData());
//	actor->GetProperty()->SetEdgeVisibility(true);
//	actor->SetMapper(mapper);
//
//	vtkNew<vtkRenderer> renderer;
//	vtkNew<vtkRenderWindow> renderWindow;
//	renderWindow->AddRenderer(renderer);
//	renderWindow->SetWindowName("CellPicking");
//
//	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//	renderWindowInteractor->Initialize();
//
//	// Set the custom stype to use for interaction.
//	vtkNew<MouseInteractorStyle> style;
//	style->SetDefaultRenderer(renderer);
//	//style->Data = triangleFilter->GetOutput();
//	style->Data = unstructuredGrid;
//
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderer->AddActor(actor);
//	renderer->ResetCamera();
//
//	renderer->SetBackground(colors->GetColor3d("PaleTurquoise").GetData());
//
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}


// An highlighted block 框选表面surface单元

//#include <vtkVersion.h>
//
//#include <vtkActor.h>
//#include <vtkAreaPicker.h>
//#include <vtkDataSetMapper.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkExtractPolyDataGeometry.h>
//#include <vtkIdFilter.h>
//#include <vtkIdTypeArray.h>
//#include <vtkInteractorStyleRubberBandPick.h>
//#include <vtkObjectFactory.h>
//#include <vtkPlanes.h>
//#include <vtkPointData.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderer.h>
//#include <vtkRendererCollection.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//#include <vtkSphereSource.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkVersion.h>
//#include <vtkVertexGlyphFilter.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkUnstructuredGridReader.h>
//#include <vtkXMLUnstructuredDataReader.h>
//#include <vtkXMLUnstructuredGridReader.h>
//#include <vtkDataSetMapper.h>
//
//
//#define VTKISRBP_ORIENT 0
//#define VTKISRBP_SELECT 1
//
//// Define interaction style
//class HighlightInteractorStyle : public vtkInteractorStyleRubberBandPick
//{
//public:
//	static HighlightInteractorStyle* New();
//	vtkTypeMacro(HighlightInteractorStyle, vtkInteractorStyleRubberBandPick);
//
//	HighlightInteractorStyle() : vtkInteractorStyleRubberBandPick()
//	{
//		this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
//		this->SelectedActor = vtkSmartPointer<vtkActor>::New();
//		this->SelectedActor->SetMapper(SelectedMapper);
//	}
//
//	virtual void OnLeftButtonUp()
//	{
//		// Forward events
//		vtkInteractorStyleRubberBandPick::OnLeftButtonUp();
//
//		if (this->CurrentMode == VTKISRBP_SELECT)
//		{
//			vtkPlanes* frustum = static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())->GetFrustum();
//
//			vtkSmartPointer<vtkExtractPolyDataGeometry> extractPolyDataGeometry =
//				vtkSmartPointer<vtkExtractPolyDataGeometry>::New();
//
//			extractPolyDataGeometry->SetInputData(this->PolyData);
//
//			extractPolyDataGeometry->SetImplicitFunction(frustum);
//			extractPolyDataGeometry->Update();
//
//			std::cout << "Extracted " << extractPolyDataGeometry->GetOutput()->GetNumberOfCells() << " cells." << std::endl;
//
//			this->SelectedMapper->SetInputData(extractPolyDataGeometry->GetOutput());
//			this->SelectedMapper->ScalarVisibilityOff();
//
//			this->SelectedActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
//			this->SelectedActor->GetProperty()->SetPointSize(5);
//
//			this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(SelectedActor);
//			this->GetInteractor()->GetRenderWindow()->Render();
//			this->HighlightProp(NULL);
//		}
//	}
//
//	void SetPolyData(vtkSmartPointer<vtkPolyData> polyData) { this->PolyData = polyData; }
//private:
//	vtkSmartPointer<vtkPolyData> PolyData;
//	vtkSmartPointer<vtkActor> SelectedActor;
//	vtkSmartPointer<vtkDataSetMapper> SelectedMapper;
//
//};
//vtkStandardNewMacro(HighlightInteractorStyle);
//
//int main(int, char* [])
//{
//	//vtkSmartPointer<vtkSphereSource> sphereSource =
//	//	vtkSmartPointer<vtkSphereSource>::New();
//	//sphereSource->Update();
//
//	vtkSmartPointer<vtkUnstructuredGridReader> reader = 
//		vtkSmartPointer<vtkUnstructuredGridReader>::New();
//	reader->SetFileName("wave1.vtk");
//	reader->Update();
//
//	auto unstructuredGrid = reader->GetOutput();
//
//	vtkSmartPointer<vtkIdFilter> idFilter =
//		vtkSmartPointer<vtkIdFilter>::New();
//	//idFilter->SetInputConnection(sphereSource->GetOutputPort());
//	idFilter->SetInputData(unstructuredGrid);
//	idFilter->Update();
//
//	// This is needed to convert the ouput of vtkIdFilter (vtkDataSet) back to vtkPolyData
//	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
//		vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//	surfaceFilter->SetInputConnection(idFilter->GetOutputPort());
//	surfaceFilter->Update();
//
//	vtkPolyData* input = surfaceFilter->GetOutput();
//
//	// Create a mapper and actor
//	vtkSmartPointer<vtkDataSetMapper> mapper =
//		vtkSmartPointer<vtkDataSetMapper>::New();
//	mapper->SetInputData(unstructuredGrid);
//	mapper->ScalarVisibilityOff();
//
//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetPointSize(5);
//
//	// Visualize
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkAreaPicker> areaPicker =
//		vtkSmartPointer<vtkAreaPicker>::New();
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetPicker(areaPicker);
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(actor);
//	//renderer->SetBackground(1,1,1); // Background color white
//
//	renderWindow->Render();
//
//	vtkSmartPointer<HighlightInteractorStyle> style =
//		vtkSmartPointer<HighlightInteractorStyle>::New();
//	style->SetPolyData(input);
//	renderWindowInteractor->SetInteractorStyle(style);
//
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}
