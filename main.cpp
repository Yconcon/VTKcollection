#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

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
//	reader->SetFileName("Test1.vtk");
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
	reader->SetFileName("zhuangpei.stl");
	reader->Update();
	input1->DeepCopy(reader->GetOutput());

	vtkNew<vtkCurvatures> curvaturesfilter;
	curvaturesfilter->SetInputConnection(reader->GetOutputPort());
	//curvaturesfilter->SetCurvatureTypeToMinimum();
	curvaturesfilter->SetCurvatureTypeToMaximum();
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
}
*/

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

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);

	renderer->SetBackground(colors->GetColor3d("Wheat").GetData());
	renderer->UseHiddenLineRemovalOn();

	vtkNew<vtkUnstructuredGridReader> reader;
	reader->SetFileName("Test1.vtk");
	reader->Update();

	//std::cout << "Loading: " << argv[1] << std::endl;
	//auto unstructuredGrid = ReadUnstructuredGrid(std::string(argv[1]));
	auto unstructuredGrid = reader->GetOutput();
	

	vtkNew<vtkLookupTable> lut1;
	lut1->SetHueRange(.667, 0);
	// Visualize
	vtkNew<vtkDataSetMapper> mapper;
	mapper->SetInputData(unstructuredGrid);
	//mapper->ScalarVisibilityOff();
	mapper->SetScalarRange(unstructuredGrid->GetScalarRange());
	mapper->SetLookupTable(lut1);
	mapper->SetColorModeToMapScalars();

	vtkNew<vtkScalarBarActor> scalarbar;
	scalarbar->SetLookupTable(mapper->GetLookupTable());
	//scalarbar->SetTitle(curvaturesfilter->GetOutput()->GetPointData()->GetScalars()->GetName());
	scalarbar->SetNumberOfLabels(5);
	renderer->AddActor2D(scalarbar);
	//vtkNew<vtkProperty> backProp;
	//backProp->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
	//backProp->SetSpecular(.6);
	//backProp->SetSpecularPower(30);

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	//actor->SetBackfaceProperty(backProp);
	//actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
	//actor->GetProperty()->SetSpecular(.3);
	//actor->GetProperty()->SetSpecularPower(30);
	//actor->GetProperty()->EdgeVisibilityOn();
	renderer->AddActor(actor);
	renderer->GetActiveCamera()->Azimuth(45);
	renderer->GetActiveCamera()->Elevation(45);
	renderer->ResetCamera();
	renderWindow->SetWindowName("ReadAllUnstructuredGridTypes");
	renderWindow->Render();
	interactor->Start();

	return EXIT_SUCCESS;
}