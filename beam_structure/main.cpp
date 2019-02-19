#define _USE_MATH_DEFINES

#include "criteria_homo.h"
#include "criteria_nonhomo.h"
#include "criteria_all.h"
#include "criteria_none.h"
#include "criteria_sheet.h"
#include "dx_surface_calculator.h"
#include "html_builder.h"
#include "matrix_fftw.h"
#include "plan_fftw.h"
#include "poly.h"
#include "povray_export.h"
#include "spacial_matrix.h"
#include "storage.h"
#include "util.h"
#include "volume_data.h"

#include "wave_matrix.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <stack>


using namespace std;
using namespace farn;

void testPresentPictures() {
  int n = 50;// Number of discrets
  double A = 0.015; // Aperture in metres
  double f = 100e6; // Frequency in Hz
  string forceProfileName("pics/sq_300_300.png");

  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(0.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  Vec3c unitForce = Vec3c(1.0, 0.0, 0.0).norm();

  try {    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitForce, FixedForce);
    MatrixFFTW forceProfile = loadFromPicture(forceProfileName, 0, 1);
    MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();

    PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    PlanFFTW backw(forceProfile, forceProfileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);
    forw.execute();

    cout << "force profile spectrum w h = " << forceProfileSpectrum.width() << " " << forceProfileSpectrum.height() << endl;
    cout << "waves dimension = " << waves.dimension() << endl;

    waves.loadFftwMatrix(forceProfileSpectrum);

#ifdef DEBUG
    waves.logState(cerr);
#endif
    
    waves.makeZShift(0.0e-3);
    Storage dat = waves.getStorage();
    Storage spacial = layerTransform(dat, forceProfile, forceProfileSpectrum, backw);
    SpacialMatrix spaceMat = getSpacialMatrixFromStorage(spacial);
    MatrixFFTW cadr(n, n);
    spaceMat.fillMatrixWithV(cadr);

    ofstream dest("picture_as_text.txt");
    dest << spaceMat << endl;

    saveAsPicture(cadr, "pics/sq_300_300_5mm_100_calculated_pic.png");
  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}

/*
  cout.precision(5);
  cout << scientific;
  
  int n = 100;
  double A = 0.001; // Aperture in metres
  double f = 60e6; // Frequency in Hz

  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(10.0 / 180.0 * M_PI);

  MaterialTensor tens = rotateTensor(tensMainAxis, rm);


  Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();

  try{
    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitForce, FixedDisplacement);
    
    MatrixFFTW forceProfile = loadFromPicture("fp_100x100.png", 0, 1);
    MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();
    PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    PlanFFTW backw(forceProfile, forceProfileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);

    forw.execute();
    waves.loadFftwMatrix(forceProfileSpectrum);

    //cout << "spectrum of force profile" << endl;
    //forceProfileSpectrum.printShortModules(cout);

    // ofstream tassLog("waves_tass.txt");
    // tassLog.precision(5);
    // tassLog << scientific;
    // waves.logState(tassLog);

    //cout << "wave matrix with amplitudes: " << endl << waves << endl << endl;

    //waves.makeZShift(0e-3);
    Storage dat = waves.getStorage();
    //ofstream dest_fs("fourier_storage.txt");
    //dat.printJson(dest_fs);
    //cout << "stored data is " << endl << dat << endl << endl;
    //cout << "fourier space storage:" << endl;
    //dat.printShortModules(cout);

    Storage spacial = layerTransform(dat, forceProfile, forceProfileSpectrum, backw);
    //cout << "spacial storage is " << endl << spacial << endl << endl;
    
    // saveAsPictures(spacial, "spacial");

    SpacialMatrix spaceMat = getSpacialMatrixFromStorage(spacial);

    //cout << "spacialMatrix = " << spaceMat << endl;
    //cout << "spacial domain storage" << endl;
    //spacial.printShortModules(cout);

    MatrixFFTW vmat(n, n);
    spaceMat.fillMatrixWithV(vmat);
    saveAsPicture(vmat, "test_pic.png");

  }catch(string msg){
    cout << "error: " << msg << endl;
  }
 */

void work() {
  int n = 600;
  double A = 0.050; // Aperture in metres
  double f = 100e6; // Frequency in Hz

  MaterialTensor tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(10.0 / 180.0 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);
  Vec3c unitDisplacement = Vec3c(1.0, 0.0, 0.0).norm();

  try{
    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitDisplacement, FixedDisplacement);
    
    string pref = "teo2_chang_10/";
    MatrixFFTW profile = loadFromPicture("fp_600_600.png", 0, 1);
    MatrixFFTW profileSpectrum = profile.cleanClone();
    PlanFFTW forw(profile, profileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    PlanFFTW backw(profile, profileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);

    forw.execute();
    waves.loadFftwMatrix(profileSpectrum);
    saveAsPicture(profileSpectrum, pref + "profile_spectrum.png");

    //waves.makeZShift(-0.0004);
    Storage dat = waves.getStorage();
    Storage spacial = layerTransform(dat, profile, profileSpectrum, backw);

    SpacialMatrix spaceMat = getSpacialMatrixFromStorage(spacial);
    Storage spaceStor(1, n, n);
    spaceMat.fillSliceWithV(0, spaceStor);
    saveAsPicture(spaceStor.sliceDW(0), pref + "from_storage.png");

    MatrixFFTW vmat(n, n);
    spaceMat.fillMatrixWithV(vmat);
    saveAsPicture(vmat, pref + "test_pic.png");

    saveAsPicture(waves.realRootsNumberMap(), pref + "rr_map.png");

  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}

void depthWork(){
  cout << "depth work" << endl;

  cout.precision(5);
  cout << scientific;

  stringstream destFlnm;
  destFlnm << "teo2_chang_10_3/"; // << setfill('0') << setw(4) << alpha / M_PI * 180;

  cout << "flnm prefix is " << destFlnm.str() << endl;  

  string suffix("_300n_105MHz_center_x_");
  
  int n = 300;// Number of discrets
  int nz = 300; // Number of z steps
  double A = 0.021; // Aperture in metres
  double f = 105e6; // Frequency in Hz
  double dz = A / n; // Shift by z in metres
  string forceProfileName("fp_300x300_center.png");

  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  //MaterialTensor tensMainAxis = makeTelluriumMaterialTensor();
  //double rho = 6210;

  //cout << "material tensor done" << endl;


  RotationMatrix rm = genChungMatrix(-4.5 / 180 * M_PI);
  //RotationMatrix rm = makeXYSlice(alpha);
  cout << "rm = " << endl << rm << endl;

  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  //double force_angle = M_PI / 4;
  Vec3c unitForce = Vec3c(1, 0, 0).norm();
  //Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();

  //cout << "chung matrix, rotate tensor and unit force are done" << endl;

  try{
    MatrixFFTW forceProfile = loadFromPicture(forceProfileName, 0, 1);
    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitForce, FixedForce);
    
    //cout << "wave matrix is done" << endl;

    
    //cout << "forceProfile is loaded" << endl;

    MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();

    //cout << "begin to make forward plan" << endl;

    PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);

    //cout << "begin to make backward plan" << endl;

    PlanFFTW backw(forceProfile, forceProfileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);

    //cout << "plans are done" << endl;

    forw.execute();

    //cout << "forward plan is executed" << endl;

    waves.loadFftwMatrix(forceProfileSpectrum);

    //cout << "fftw matrix is loaded" << endl;

    //cout << "try to create volume storage h = " << nz + 1 << " d = " << n << " w = " << n << endl;
    
    Storage volumeStorage(nz + 1, n, n);

    //cout << "volume storage created" << endl;

    for(int t = 0; t <= nz; ++t){
      cout << "t = " << t << endl;

      Storage storageFourier = waves.getStorage();
      Storage storageSpacial = layerTransform(storageFourier, forceProfile, forceProfileSpectrum, backw);
      SpacialMatrix spaceMat = getSpacialMatrixFromStorage(storageSpacial);

      spaceMat.fillSliceWithV(t, volumeStorage);

      waves.makeZShift(-dz);
    }

    volumeStorage.save(destFlnm.str() + suffix + "volumeData.txt");
    
    MatrixFFTW sliceZX = volumeStorage.sliceHW(n / 2);
    MatrixFFTW sliceZY = volumeStorage.sliceHD(n / 2);

    saveAsPicture(sliceZX, destFlnm.str() + suffix + "zx.png");
    saveAsPicture(sliceZY, destFlnm.str() + suffix + "zy.png");
    

    //volumeStorage.saveForDx("15mhz_zpol_te_30x11_volume_data");

    //cout << "volume storage is " << endl << volumeStorage << endl;
  } catch(string msg) {
    cout << "error: " << msg << endl;
  }
}

void testPol(){
  cerr << "testPol" << endl;
  double sq = sqrt(2.0);

  ComplexVec v1 = list_of (0) (0) (1);
  ComplexVec v2 = list_of (1/sq) (1/sq) (0);
  ComplexVec v3 = list_of (0.5) (0.5) (1/sq);
  ComplexVec zero = list_of (0) (0) (0);

  cout << vecAbs(v1) << " " << vecAbs(v2) << " " << vecAbs(v3) << endl;

  cout << calcPol(matrixOfRows(v1,v2,v3)) << endl << endl;
  cout << calcPol(matrixOfRows(v1,v1,v3)) << endl << endl;
  cout << calcPol(matrixOfRows(v1,v2,v2)) << endl << endl;
  cout << calcPol(matrixOfRows(v2,v3,v3)) << endl << endl;
  cout << calcPol(matrixOfRows(v1,v2,v1)) << endl << endl;
  cout << calcPol(matrixOfRows(v2,v2,v3)) << endl << endl;
  cout << calcPol(matrixOfRows(v3,v2,v3)) << endl << endl;

  cout << calcPol(matrixOfRows( 
			       (list_of (-7.95e10) (0.0) (0.0) ), 
			       (list_of (0.0) (-7.95e10) (0.0) ), 
			       (list_of (0.0) (0.0) (0.0) ))) << endl << endl;

  
}

void testPol1() {
  CD mat[3][3] = {
    {CD(-1102.84,0), CD(-4.0958e-13,0), CD(1.24549e-13,0)}, 
    {CD(-5.84609e-15,0), CD(14091.7,0), CD(71.8375,0)}, 
    {CD(1.37198e-13,0), CD(71.8375,0), CD(0.366217,0)}};
  
  Mat3 tm;
  for(int p = 0; p < 3; p++) {
    for(int q = 0; q < 3; q++) {
      tm[p][q] = mat[p][q];
    }
  }

  Vec3c pol = tm.calcPol();
  cout << "tm = " << endl;

  for(int p = 0; p < 3; p++) {
    for(int q = 0; q < 3; q++) {
      cout << tm[p][q] << " ";
    }
    cout << endl;
  }

  cout << "pol = " << pol << endl;
  cout << "tm * pol = " << (tm * pol) << endl;
}

void testCompositeWave(){
  cerr << "begin of test composite wave" << endl;
  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(10.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);
  cerr << "rotate tensor" << endl;

  Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();
  cerr << "unitForce = " << unitForce << endl;

  try{
    CompositeWave cv(0.0, -0.00044, tens, rho, 2 * M_PI * 1e8, unitForce, FixedForce);
    cerr << "cv done" << endl;

    cout << endl;
    cout << cv << endl;
  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}

void testCompositeWave1(){
  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(9.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  Vec3c unitForce = Vec3c(0.2, 0.3, 0.4).norm();
  cerr << "unitForce = " << unitForce << endl;

  try{
    CompositeWave cv(0.0, 0.0, tens, rho, 2 * M_PI * 1e8, unitForce, FixedForce);
    cout << endl;
    cout << cv << endl;
  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}


void
testStorage(){
  int h = 3, d = 4, w = 5;

  Storage dat(h, d, w);

  for(int t = 0; t < h; ++t)
    for(int p = 0; p < d; ++p)
      for(int q = 0; q < w; ++q)
	dat(t, p, q) = (t + 1) * 100 + (p + 1) * 10 + (q + 1);

  cout << "storage: " << endl << dat << endl << endl;
    
}

SpacialMatrix 
makeTestSpacialMatrix(int n){
  SpacialMatrix sm(n, n);

  for(int p = 0; p < n; ++p) {
    for(int q = 0; q < n; ++q) {
      sm(p,q).setTestValue(p, q);
    }
  }

  return sm;
}

void
testStorageSpacialMatrix(){
  int n = 10;

  SpacialMatrix sm = makeTestSpacialMatrix(n);

  cout << "spacial matrix " << endl << sm << endl << endl;

  Storage stor(1, n, n);
  sm.fillSliceWithV(0, stor);

  cout << "storage: " << endl << stor << endl;
}

void
testRootSelectionAlgo() {
  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(9.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  Vec3c unitForce = Vec3c(0.2, 0.3, 0.4).norm();
  cerr << "unitForce = " << unitForce << endl;

  try{
    double smin = -0.002;
    double smax =  0.002;

    int n = 1000;
    double ds = (smax - smin) / n;
    
    for(double sx = smin; sx < smax; sx += ds) {
      for(double sy = smin; sy < smax; sy += ds) {
	CompositeWave cv(sx, sy, tens, rho, 2 * M_PI * 1e8, unitForce, FixedForce);
      }
    }
  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}

void testStorageJson() {
  int h = 2;
  int d = 3;
  int w = 4;

  Storage stor(h, d, w);
  for(int p = 0; p < h; ++p) {
    for(int q = 0; q < d; ++q) {
      for(int r = 0; r < w; ++r) {
	stor(p, q, r) = 1000 + 100 * p + 10 * q + r;
      }
    }
  }

  stor.printJson(cout);
}

void testSelection() {
  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;
  cout.precision(3);
  try{
    int n = 10;
    double A = 0.02;
    double f = 20e6;
    Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();
    WaveMatrix waves = create_wave_matrix(n, A, f, tensMainAxis, rho, unitForce, FixedForce);
    waves.setOnes();

    WaveMatrix homo = waves;
    WaveMatrix nonHomo = waves;

    homo.filter(CriteriaNonHomo());
    nonHomo.filter(CriteriaHomo());

    for(int number = 0; number < 15; number ++) {
      stringstream name1, name2;
      name1 << "pic_homo_" << number << ".png";
      name2 << "pic_nonhomo_" << number << ".png";
    
      MatrixFFTW pic1 = homo.getStorage().sliceDW(number);
      MatrixFFTW pic2 = nonHomo.getStorage().sliceDW(number);
      
      //cout << "homo " << number << endl << pic1 << endl;
      //cout << "nonhomo " << number << endl << pic2 << endl;

      //sav eAsPicture(pic1, name1.str().c_str());
      //saveAsPicture(pic2, name2.str().c_str());
    }
  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}

void testCompositeWaveZ() {
  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;
  //double A = 0.02;
  //double f = 20e6;
  Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();
  try {
    CompositeWave cw(0, 0, tensMainAxis, rho, 20e8, unitForce, FixedForce);
  } catch(string msg) {
    cout << "exception: " << msg << endl;
  }
}

void freqSeq() {
  cerr << "let the game begin!" << endl;
  double rho = 5.96e3;//material density
  double A = 0.001; //Aperture in metres
  int n = 300;// Number of discrets
  int nz = 300; // Number of z steps
  double dz = A / n;
  Vec3c unitDisplacement = Vec3c(0.0, 0.0, 1.0).norm(); //force direction
  string pref = "res/";

  string forceProfileName("pics/sq_300_300_5x5.png");
  MatrixFFTW forceProfile = loadFromPicture(forceProfileName, 0, 1);
  MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();
  PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
  PlanFFTW backw(forceProfile, forceProfileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);
  forw.execute();
  MatrixFFTW saveForceProfileSpectrum = forceProfileSpectrum;

  MaterialTensor tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  RotationMatrix rm = genChungMatrix(10.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  Storage volumeStorage(nz + 1, n, n);
  Storage volumeStorageHomo(nz + 1, n, n);
  Storage volumeStorageNonHomo(nz + 1, n, n);

  HtmlBuilder table("polz_disp.html");
  table.beginRow();
  table.setHeaderItem("fourier");
  table.setHeaderItem("roots map");
  table.setHeaderItem("eliminate map");

  table.setHeaderItem("b zero z");
  table.setHeaderItem("b half z");
  table.setHeaderItem("b zx");
  table.setHeaderItem("b zy");

  table.setHeaderItem("h zero z");
  table.setHeaderItem("h half z");
  table.setHeaderItem("h zx");
  table.setHeaderItem("h zy");

  table.setHeaderItem("n half z");
  table.setHeaderItem("n zero z");
  table.setHeaderItem("n zx");
  table.setHeaderItem("n zy");
  table.endRow();

  cerr << "ready to cycle" << endl;

  for(double fshow = 20; fshow <= 220; fshow += 20) {
    cerr << " f = " << fshow << endl;
    double f = fshow * 1e6;
    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitDisplacement, FixedDisplacement);
    waves.loadFftwMatrix(saveForceProfileSpectrum);

    for(int t = 0; t <= nz; ++t) {
      cout << "\t t = " << t << endl;
      WaveMatrix waves_homo = waves;
      WaveMatrix waves_non_homo = waves;
      waves_homo.filter(CriteriaHomo());
      waves_non_homo.filter(CriteriaNonHomo());

      Storage storageFourier = waves.getStorage();
      Storage homoFourier = waves_homo.getStorage();
      Storage nonHomoFourier = waves_non_homo.getStorage();

      Storage storageSpacial = layerTransform(storageFourier, forceProfile, forceProfileSpectrum, backw);
      Storage storageSpacialHomo = layerTransform(homoFourier, forceProfile, forceProfileSpectrum, backw);
      Storage storageSpacialNonHomo = layerTransform(nonHomoFourier, forceProfile, forceProfileSpectrum, backw);

      SpacialMatrix spaceMat = getSpacialMatrixFromStorage(storageSpacial);
      SpacialMatrix spaceMatHomo = getSpacialMatrixFromStorage(storageSpacialHomo);
      SpacialMatrix spaceMatNonHomo = getSpacialMatrixFromStorage(storageSpacialNonHomo);

      spaceMat.fillSliceWithV(t, volumeStorage);
      spaceMatHomo.fillSliceWithV(t, volumeStorageHomo);
      spaceMatNonHomo.fillSliceWithV(t, volumeStorageNonHomo);

      waves.makeZShift(-dz);
    }

    MatrixFFTW sliceZX = volumeStorage.sliceHW(n / 2);
    MatrixFFTW sliceZY = volumeStorage.sliceHD(n / 2);
    MatrixFFTW sliceHalfZ = volumeStorage.sliceDW(nz / 2);
    MatrixFFTW sliceZeroZ = volumeStorage.sliceDW(0);

    MatrixFFTW sliceZXHomo = volumeStorageHomo.sliceHW(n / 2);
    MatrixFFTW sliceZYHomo = volumeStorageHomo.sliceHD(n / 2);
    MatrixFFTW sliceHalfZHomo = volumeStorageHomo.sliceDW(nz / 2);
    MatrixFFTW sliceZeroZHomo = volumeStorageHomo.sliceDW(0);

    MatrixFFTW sliceZXNonHomo = volumeStorageNonHomo.sliceHW(n / 2);
    MatrixFFTW sliceZYNonHomo = volumeStorageNonHomo.sliceHD(n / 2);
    MatrixFFTW sliceHalfZNonHomo = volumeStorageNonHomo.sliceDW(nz / 2);
    MatrixFFTW sliceZeroZNonHomo = volumeStorageNonHomo.sliceDW(0);

    MatrixFFTW realRootsNumberMap = waves.realRootsNumberMap();

    stringstream flnmFourier, flnmMap, flnmEliminateMap;

    stringstream flnmPlaneZeroZ;
    stringstream flnmPlaneZX;
    stringstream flnmPlaneZY;
    stringstream flnmPlaneHalfZ;

    stringstream flnmHomoPlaneZeroZ;
    stringstream flnmHomoPlaneZX;
    stringstream flnmHomoPlaneZY;
    stringstream flnmHomoPlaneHalfZ;

    stringstream flnmNonHomoPlaneZeroZ;
    stringstream flnmNonHomoPlaneZX;
    stringstream flnmNonHomoPlaneZY;
    stringstream flnmNonHomoPlaneHalfZ;
   
    flnmFourier << pref << fshow << "_fourier.png";
    flnmMap << pref << fshow << "_map.png";
    flnmEliminateMap << pref << fshow << "_eliminate_map.png";

    flnmPlaneZeroZ << pref << fshow << "_b_zero_z.png";
    flnmPlaneZX << pref << fshow << "_b_zx.png";
    flnmPlaneZY << pref << fshow << "_b_zy.png";
    flnmPlaneHalfZ << pref << fshow << "_b_half_z.png";

    flnmHomoPlaneZeroZ << pref << fshow << "_h_zero_z.png";
    flnmHomoPlaneZX << pref << fshow << "_h_zx.png";
    flnmHomoPlaneZY << pref << fshow << "_h_zy.png";
    flnmHomoPlaneHalfZ << pref << fshow << "_h_half_z.png";

    flnmNonHomoPlaneZeroZ << pref << fshow << "_n_zero_z.png";
    flnmNonHomoPlaneZX << pref << fshow << "_n_zx.png";
    flnmNonHomoPlaneZY << pref << fshow << "_n_zy.png";
    flnmNonHomoPlaneHalfZ << pref << fshow << "_n_half_z.png";

    saveAsPicture(saveForceProfileSpectrum, flnmFourier.str());
    saveAsPicture(realRootsNumberMap,   flnmMap.str());
    saveAsPicture(waves.eliminateMap(), flnmEliminateMap.str());
    
    saveAsPicture(sliceZeroZ, flnmPlaneZeroZ.str());
    saveAsPicture(sliceZX, flnmPlaneZX.str());
    saveAsPicture(sliceZY, flnmPlaneZY.str());
    saveAsPicture(sliceHalfZ, flnmPlaneHalfZ.str());

    saveAsPicture(sliceZeroZHomo, flnmHomoPlaneZeroZ.str());
    saveAsPicture(sliceZXHomo, flnmHomoPlaneZX.str());
    saveAsPicture(sliceZYHomo, flnmHomoPlaneZY.str());
    saveAsPicture(sliceHalfZHomo, flnmHomoPlaneHalfZ.str());

    saveAsPicture(sliceZeroZNonHomo, flnmNonHomoPlaneZeroZ.str());
    saveAsPicture(sliceZXNonHomo, flnmNonHomoPlaneZX.str());
    saveAsPicture(sliceZYNonHomo, flnmNonHomoPlaneZY.str());
    saveAsPicture(sliceHalfZNonHomo, flnmNonHomoPlaneHalfZ.str());

    table.beginRow();
    table.setPicture(flnmFourier.str());
    table.setPicture(flnmMap.str());
    table.setPicture(flnmEliminateMap.str());

    table.setPicture(flnmPlaneZeroZ.str());
    table.setPicture(flnmPlaneHalfZ.str());
    table.setPicture(flnmPlaneZX.str());
    table.setPicture(flnmPlaneZY.str());

    table.setPicture(flnmHomoPlaneZeroZ.str());
    table.setPicture(flnmHomoPlaneHalfZ.str());
    table.setPicture(flnmHomoPlaneZX.str());
    table.setPicture(flnmHomoPlaneZY.str());

    table.setPicture(flnmNonHomoPlaneZeroZ.str());
    table.setPicture(flnmNonHomoPlaneHalfZ.str());
    table.setPicture(flnmNonHomoPlaneZX.str());
    table.setPicture(flnmNonHomoPlaneZY.str());
    table.endRow();
  }
}

void testFreqSeq() {
  try {
    freqSeq();
  } catch(string msg) {
    cerr << "error: " << msg << endl;
  }
}

void testPureFourier() {
  string forceProfileName("pics/sq_300_300.png");
  MatrixFFTW forceProfile = loadFromPicture(forceProfileName, 0, 1);
  MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();
  PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
  forw.execute();
  saveAsPicture(forceProfile, "test_source.png");
  saveAsPicture(forceProfileSpectrum, "test_destination.png");
}

void testHtmlBuilder() {
  cerr << "test html builder" << endl;


  HtmlBuilder hb("test.html");
  hb.beginRow();
  hb.setHeaderItem("a");
  hb.setHeaderItem("b");
  hb.setHeaderItem("c");
  hb.endRow();

  hb.beginRow();
  hb.setPicture("11.png");
  hb.setPicture("12.png");
  hb.setPicture("13.png");
  hb.endRow();

  hb.beginRow();
  hb.setPicture("21.png");
  hb.setPicture("22.png");
  hb.setPicture("23.png");
  hb.endRow();

}

void testTransmission() {
  cout << "test transmission" << endl;
  int n = 500;
  double A = 0.005; // Aperture in metres
  double f = 24e6; // Frequency in Hz

  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(10.0 / 180.0 * M_PI);

  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();

  //0.0000643852, 0.0000643855, 0.0005258366, 0.0005258367
  //0.0, 0.002, 0.00, 0.002,
  WaveMatrix wm = create_framed_wave_matrix(n, 0.0, 0.002, 0.00, 0.002,  A, f, tens, rho,  unitForce, FixedForce);
  MatrixFFTW mat = wm.transmission();
  saveAsPicture(mat, "transmission.png");

  MatrixFFTW transmissionSpectrum = mat.cleanClone();
  PlanFFTW forw(mat, transmissionSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
  forw.execute();

  saveAsPicture(transmissionSpectrum, "transmission_spectrum.png");
  ofstream dest("transmission_matrix.txt");
  dest << mat;
  cout << "the end" << endl;
}

void testNan() {
  double rho = 5.96e3;//material density
  double A = 0.001; //Aperture in metres
  double f = 20e6; //frequency Hz
  int n = 2000;// Number of discrets
  Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm(); //force direction

  MaterialTensor tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  RotationMatrix rm = genChungMatrix(-10.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  WaveMatrix wm = create_wave_matrix(n, A, f, tens, rho, unitForce, FixedForce);
}

void testDisplacement() {
  double rho = 5.96e3;//material density
  double f = 20e6; //frequency Hz
  Vec3c unitDisplacement = Vec3c(0.0, 1.0, 0.0).norm(); //force direction

  MaterialTensor tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  RotationMatrix rm = genChungMatrix(-10.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  for(double sx = -0.006; sx <= 0.006; sx += 0.00001) {
    for(double sy = -0.006; sy <= 0.006; sy += 0.00001) {
      CompositeWave cw(sx, sy, tens, rho, f, unitDisplacement, FixedDisplacement);
      cw.setAmplitudeOne();

      Storage stor(15, 1, 1);
      cw.incrementStorage(stor, 0, 0);
      Vec3c res(stor(0, 0, 0), stor(1, 0, 0), stor(2, 0, 0));
      if ((res - unitDisplacement).abs() > 1e-10) {
	cout << "break at s1 = " << sx << ", s2 = " << sy << endl;
      }
    }
  }
  cout << "done" << endl;
}

mutex splitModesMutex;

void threadSplitModes(stack<int>* tasks, double dz, WaveMatrix *pwm1, WaveMatrix *pwm2, WaveMatrix *pwm3, WaveMatrix *pwm4,
		      Storage* pVolumeStorage, Storage* pVolumeStorage1, Storage* pVolumeStorage2, Storage* pVolumeStorage3 ) {
  
  MatrixFFTW source(pVolumeStorage->width(), pVolumeStorage->depth());
  MatrixFFTW dest = source.cleanClone();

  splitModesMutex.lock();
  PlanFFTW backw(source, dest, FFTW_BACKWARD, FFTW_ESTIMATE);
  splitModesMutex.unlock();

  while(tasks->size()) {
    int task;
    splitModesMutex.lock();
    task = tasks->top();
    tasks->pop();
    splitModesMutex.unlock();

    cerr << task << endl;

    WaveMatrix waves = *pwm1;
    WaveMatrix waves_mode0 = *pwm2;
    WaveMatrix waves_mode1 = *pwm3;
    WaveMatrix waves_mode2 = *pwm4;

    waves.makeZShift(-dz*task);
    waves_mode0.makeZShift(-dz*task);
    waves_mode1.makeZShift(-dz*task);
    waves_mode2.makeZShift(-dz*task);

    Storage storageFourier = waves.getStorage();
    Storage storageSpacial = layerTransform(storageFourier, source, dest, backw);
    SpacialMatrix spaceMat = getSpacialMatrixFromStorage(storageSpacial);
    spaceMat.fillSliceWithV(task, *pVolumeStorage);

    Storage storageFourier0 = waves_mode0.getStorage();
    Storage storageSpacial0 = layerTransform(storageFourier0, source, dest, backw);
    SpacialMatrix spaceMat0 = getSpacialMatrixFromStorage(storageSpacial0);
    spaceMat0.fillSliceWithV(task, *pVolumeStorage1);

    Storage storageFourier1 = waves_mode1.getStorage();
    Storage storageSpacial1 = layerTransform(storageFourier1, source, dest, backw);
    SpacialMatrix spaceMat1 = getSpacialMatrixFromStorage(storageSpacial1);
    spaceMat1.fillSliceWithV(task, *pVolumeStorage2);

    Storage storageFourier2 = waves_mode2.getStorage();
    Storage storageSpacial2 = layerTransform(storageFourier2, source, dest, backw);
    SpacialMatrix spaceMat2 = getSpacialMatrixFromStorage(storageSpacial2);
    spaceMat2.fillSliceWithV(task, *pVolumeStorage3);
  }
}

void testSplitModes() {
  cerr << "test split modes" << endl;
  double rho = 5.96e3;//material density
  double A = 0.010; //Aperture in metres
  int n = 300;// Number of discrets
  int nz = 150; // Number of z steps
  double dz = A / n;
  double f = 130e6;
  Vec3c unitDisplacement = Vec3c(0.0, 0.0, 1.0).norm(); //force direction
  
  MaterialTensor tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  RotationMatrix rm = nearXMatrix(4.0 / 180 * M_PI);
  string pref = "split_modes00/4_deg_polx_";
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitDisplacement, FixedForce);

  string forceProfileName("fp_300x300_center.png"/*"pics/sq_300_300.png"*//*"fp_100x100.png"*/);
  MatrixFFTW forceProfile = loadFromPicture(forceProfileName, 0, 1);
  MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();
  PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
  PlanFFTW backw(forceProfile, forceProfileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);
  forw.execute();
  waves.loadFftwMatrix(forceProfileSpectrum);

  WaveMatrix waves_mode0 = waves;
  WaveMatrix waves_mode1 = waves;
  WaveMatrix waves_mode2 = waves;

  waves_mode0.filter(CriteriaSheet(0));
  waves_mode1.filter(CriteriaSheet(1));
  waves_mode2.filter(CriteriaSheet(2));

  WaveMatrix wm0 = waves_mode0;
  WaveMatrix wm1 = waves_mode1;
  WaveMatrix wm2 = waves_mode2;

  Storage volumeStorage(nz + 1, n, n);
  Storage volumeStorage0(nz + 1, n, n);
  Storage volumeStorage1(nz + 1, n, n);
  Storage volumeStorage2(nz + 1, n, n);

  stack<int> tasks;
  for(int t = 0; t <= nz; ++t) {
    tasks.push(t);
  }

  thread threads[4];
  for(int t = 0; t < 4; ++t) {
    threads[t] = thread(threadSplitModes, &tasks, dz, &waves, 
			&waves_mode0, &waves_mode1, &waves_mode2, 
			&volumeStorage, &volumeStorage0, &volumeStorage1, &volumeStorage2);
  }

  for(int t = 0; t < 4; ++t) {
    threads[t].join();
  }


  saveAsPicture(volumeStorage.sliceHW(n / 2),  pref + "slice_zx.png");
  saveAsPicture(volumeStorage.sliceHD(n / 2),  pref + "slice_zy.png");
  saveAsPicture(volumeStorage.sliceDW(0.05 * nz),  pref + "slice_xy_near.png");
  saveAsPicture(volumeStorage.sliceDW( 0.5 * nz),  pref + "slice_xy_far.png");

  saveAsPicture(volumeStorage0.sliceHW(n / 2),  pref + "0_slice_zx.png");
  saveAsPicture(volumeStorage0.sliceHD(n / 2),  pref + "0_slice_zy.png");
  saveAsPicture(volumeStorage0.sliceDW(0.05 * nz),  pref + "0_slice_xy_near.png");
  saveAsPicture(volumeStorage0.sliceDW( 0.5 * nz),  pref + "0_slice_xy_far.png");

  saveAsPicture(volumeStorage1.sliceHW(n / 2),  pref + "1_slice_zx.png");
  saveAsPicture(volumeStorage1.sliceHD(n / 2),  pref + "1_slice_zy.png");
  saveAsPicture(volumeStorage1.sliceDW(0.05 * nz),  pref + "1_slice_xy_near.png");
  saveAsPicture(volumeStorage1.sliceDW( 0.5 * nz),  pref + "1_slice_xy_far.png");

  saveAsPicture(volumeStorage2.sliceHW(n / 2),  pref + "2_slice_zx.png");
  saveAsPicture(volumeStorage2.sliceHD(n / 2),  pref + "2_slice_zy.png");
  saveAsPicture(volumeStorage2.sliceDW(0.05 * nz),  pref + "2_slice_xy_near.png");
  saveAsPicture(volumeStorage2.sliceDW( 0.5 * nz),  pref + "2_slice_xy_far.png");
}

void wrapSplitModes() {
  try {
    testSplitModes();
  } catch(string msg) {
    cout << "error: " << msg << endl;
  }
}

void workTellurium() {
  MaterialTensor mt = makeTelluriumMaterialTensor();
  double rho = 6210;

  Vec3 a=Vec3(cos(10*M_PI/180), sin(10*M_PI/180), 0).norm(), b(0, 0, 1);
  int n = 100;
  double dphi = 2 * M_PI / n;
  double dz = 2.0 / n;

  for(double z = -1; z <= 1; z += dz) {
    for(double phi = 0; phi < 2 * M_PI; phi += dphi) {
      double r = sqrt(1 - z * z);
      Vec3 vc(r * cos(phi), r * sin(phi), z);
      RotationMatrix cm = makeChristoffelMatrix(mt, vc);
      Poly pm = makeCharacterPoly(cm);
      Poly::RootVec ar = pm.all_roots();
      vector<double> vels(ar.size());
      transform(ar.begin(), ar.end(), vels.begin(), [rho](complex<double>& x) -> double {return sqrt(abs(x)/rho);});
      sort(vels.begin(), vels.end());
      for(auto it = ar.begin(); it != ar.end(); ++it) {
	double vel = sqrt(abs(*it) / rho);
	double xv = vc[0] / vel;
	double yv = vc[1] / vel;
	double zv = vc[2] / vel;

	cout << xv << " " << yv << " " << zv << endl;
      }
    }
  }
}

void workMainTelluriumDirections() {
  MaterialTensor mt = makeTelluriumMaterialTensor();
  double rho = 6210;

  Vec3 vc=Vec3(1, 1, 0).norm();
  RotationMatrix cm = makeChristoffelMatrix(mt, vc);
  Poly pm = makeCharacterPoly(cm);
  Poly::RootVec ar = pm.all_roots();
  vector<double> vels(ar.size());
  transform(ar.begin(), ar.end(), vels.begin(), [rho](complex<double>& x) -> double {return sqrt(abs(x)/rho);});
  sort(vels.begin(), vels.end());
  for(auto it = vels.begin(); it != vels.end(); ++it) {
    cout << "\t" << (*it) << endl;
  }
}

void dxSlowness() {
  MaterialTensor mt = makeTelluriumMaterialTensor();
  double rho = 6210;
  //MaterialTensor mt = makeParatelluriteMaterialTensor();
  //double rho = 5970;

  DxSurfaceCalculator calc(400, mt, rho);
  calc.save("te_pol_3", 2);
}

void workColorSheets() {
  cerr << "color sheets are in process" << endl;
  int n = 300;// Number of discrets
  //int nz = 150; // Number of z steps
  double A = 0.08; // Aperture in metres
  double f = 20.0e6; // Frequency in Hz
  //double dz = A / n; // Shift by z in metres
  string forceProfileName("pics/sq_300_300_30x11.png");
  //string forceProfileName("pics/sq_300_300_5x5.png");
  //string forceProfileName("fp_10x10.png");

  //MaterialTensor tensMainAxis = makeParatelluriteMaterialTensor();
  //double rho = 5970;

  MaterialTensor  tensMainAxis = makeTelluriumMaterialTensor();
  double rho = 6210;

  double cutAngle = 5;
  RotationMatrix rm = makeAxisRotation(M_PI/2, 1); //genChungMatrix(10.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(rotateTensor(tensMainAxis, makeAxisRotation(-cutAngle/180*M_PI, 2)), makeAxisRotation(90.0/180*M_PI, 1));
  Vec3c unitForce = Vec3c(0.0, 0.0, 1.0).norm();

  try{
    MatrixFFTW forceProfile = loadFromPicture(forceProfileName, 0, 1);
    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitForce, FixedDisplacement);

    MatrixFFTW forceProfileSpectrum = forceProfile.cleanClone();
    PlanFFTW forw(forceProfile, forceProfileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    PlanFFTW backw(forceProfile, forceProfileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);
    forw.execute();
    waves.loadFftwMatrix(forceProfileSpectrum);
    ofstream destSheets("sheets_te_along_x+5_pol_z.dx");
    waves.dxExport(destSheets);
  } catch(string msg) {
    cout << "error: " << msg << endl;
  }
}

double distFunc(const MaterialTensor& tens, double rho, const Vec3& a, const Vec3& b, double phi, int ind) {
    Vec3 vc = a * cos(phi) + b * sin(phi);
    RotationMatrix cm = makeChristoffelMatrix(tens, vc);
    Poly pm = makeCharacterPoly(cm);
    Poly::RootVec ar = pm.all_roots();
    vector<double> vels(ar.size());
    transform(ar.begin(), ar.end(), vels.begin(), [rho](complex<double>& x) -> double {return sqrt(abs(x)/rho);});
    sort(vels.begin(), vels.end());
    return 1/vels[ind];
}

void 
workWalkOff() {
  MaterialTensor tensMainAxis = makeParatelluriteMaterialTensor();
  double rho = 5970;
  ofstream dest("walkoff_2.txt");
  int ind = 2;

  Vec3 a=Vec3(1, 0, 0).norm();
  Vec3 b=Vec3(0, 1, 0).norm();

  int n = 10000;
  double dphi = (M_PI/2) / n;
  double val1 = 0;
  double val2 = 0;
  bool flag = true;

  for(double phi = 0; phi < M_PI/2; phi += dphi) {
    val1 = val2;
    val2 = distFunc(tensMainAxis, rho, a, b, phi, ind);
    if (flag) {
      flag = false;
    } else {
      double deriv = (val2 - val1) / dphi;
      dest << phi/M_PI*180 << " " << atan2(deriv, (val1 + val2) / 2)/M_PI*180 << endl;
    }
  }
  
}

void testRoots() {
  Poly sam1((Poly::CoeffVec const&)(list_of(-2.16e11)(0)(5.724e8)(0)(-379215)(0)(74.4385)));
  vector<complex<double> > allRoots = sam1.all_roots();
  cout << "roots are " << endl;
  cout << allRoots << endl;
}

void testPolynomial() {
  MaterialTensor mt = makeParatelluriteMaterialTensor();
  double rho = 5980;


  double s1 = 0.001;
  double s2 = 0.0002;

  PolyMatrix pm = makePolyMatrix(mt, rho, s1,  s2);

  cout << "s1 = " << s1 << " s2 = " << s2 << endl;

  cout << "poly matrix = " << endl << pm << endl;

  Poly p = det(pm);

  cout << "polynome = " << p << endl;

  Poly::RootVec ar = p.all_roots();

  cout << "roots are " << ar << endl;

  Vec3 sf1(s1, s2, real(ar[0]));
  Vec3 sf2(s1, s2, real(ar[1]));
  Vec3 sf3(s1, s2, real(ar[2]));
  Vec3 sf4(s1, s2, real(ar[3]));
  Vec3 sf5(s1, s2, real(ar[4]));
  Vec3 sf6(s1, s2, real(ar[5]));

  cout << "v1 = " << 1/sf1.abs() << endl;
  cout << "v2 = " << 1/sf2.abs() << endl;
  cout << "v3 = " << 1/sf3.abs() << endl;
  cout << "v4 = " << 1/sf4.abs() << endl;
  cout << "v5 = " << 1/sf5.abs() << endl;
  cout << "v6 = " << 1/sf6.abs() << endl;

}

void testChang() {
  MaterialTensor  tensMainAxis = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 2.65e10, 6.6e10);
  double rho = 5.96e3;

  RotationMatrix rm = genChungMatrix(0.0 / 180 * M_PI);
  MaterialTensor tens = rotateTensor(tensMainAxis, rm);

  ofstream dest("slow_surf.txt");

  double a = 0.002;
  int n = 100;

  for(double sy = -a; sy < a; sy += (a/n)) {
    for(double sx = -a; sx < a; sx += (a/n)) {
	PolyMatrix poly_mat = makePolyMatrix(tens, rho, sx, sy);
	Poly pol = det(poly_mat);
	Poly::RootVec roots = pol.all_roots();

	bool flag = true;
	double val = 0;
	for(size_t ind = 0; ind < roots.size(); ++ind) {
	  if (abs(imag(roots[ind])) < 1e-6) {
	    if (flag) {
	      flag = false;
	      val = real(roots[ind]);
	    } else {
	      val = max(val, real(roots[ind]));
	    }
	  }
	}

	dest << val << " ";
    }
    dest << endl;
  }
}

double index2value(int r, int p, int q) {
  return sin(sqrt(r*r + p*p + q*q));
}

void testStorageSaveLoad() {
  int h = 20, d = 7, w = 12;

  Storage a(h, d, w), b;
  
  for(int r = 0; r < h; r++) {
    for(int p = 0; p < d; p++) {
      for(int q = 0; q < w; q++) {
	a(r, p , q) = index2value(r, p, q);
      }
    }
  }

  a.save("test_storage.txt");

  b.load("test_storage.txt");

  if (a==b) {
    cout << "it is OK with the save/load functions" << endl;
  } else {
    for(int r = 0; r < h; r++) {
      for(int p = 0; p < d; p++) {
	for(int q = 0; q < w; q++) {
	  if (abs(a(r, p, q) - b(r, p, q)) > 1e-8) {
	      cerr << "mismatch in " << r << " " << p << " " << q << endl;
	      cerr << "\t" << a(r, p, q) << " != " << b(r, p, q) << endl;
	          
	      cout << "something wrong with the serialization" << endl;
	  }
	}
      }
    }
  }
}

void testBrezen() {
  int r = 0, p = 0, q = 0;
  Vec3 dir(0.3, 0.4, 0.5);

  PovrayExport povBrezen("povbrezen.pov");

  while(r < 300) {
    bool firstTime = true;
    double deviat = 0;
    int rr = 0, pp = 0, qq = 0; 

    for(int dp = -1; dp < 2; dp++) {
      for(int dq = -1; dq < 2; dq++) {
	for(int dr = 0; dr < 2; dr++) {
	  if (dp != 0 && dq != 0 && dr != 0) {
	    int nr = r + dr, np = p + dp, nq = q + dq;
	    double val = farn::brezenDeviat(Vec3(nr, np, nq), dir);
	    if (firstTime) {
	      firstTime = false;
	      deviat = val;
	      rr = nr;
	      pp = np; 
	      qq = nq;
	    } else {
	      if (val < deviat) {
		deviat = val;
		pp = np;
		qq = nq;
		rr = nr;
	      }
	    }
	  }
	}
      }
    }
    
    povBrezen.sphere(r, p, q);

    cerr << p << " " <<  q << " " << r << endl;
    r = rr;
    p = pp; 
    q = qq;
  }
}

void testAverageStorage() {
  Storage volume;
  volume.load("teo2_chang_10_3/_300n_105MHz_center_x_volumeData.txt");

  cerr << "size = " << volume.height() << " " << volume.depth() << " " << volume.width() << endl;

  double wa_angle = -38 * M_PI / 180;
  MatrixFFTW mat = volume.accumulateAlongDirection(Vec3(0, sin(wa_angle), cos(wa_angle)));
  saveAsPicture(mat, "projection.png");
}


void testVolumeVisualization() {
  Mat3 cs(
	  Vec3c(0, 0, 1), 
	  Vec3c(1, 0, 0), 
	  Vec3c(0, 1, 0));

  VolumeData svd = newBarZVolumeData(100, cs);
  drawVolumePovray("spherical.pov", svd, 0.97);
}

void workTelluriumBeamStructure() {
  int n = 600;
  double A = 0.013; // Aperture in metres
  double f = 20e6;  // Frequency in Hz
  double dz = A / n;

  //cerr << "transducer dx in pix = " << n * (0.006/A) << endl;
  //cerr << "transducer dy in pix = " << n * (0.006/A) << endl;
  //return;

  MaterialTensor tensMainAxis = makeTelluriumMaterialTensor();
  double rho = 6210;

  RotationMatrix zxTransdMatr(combineInMatrix(Vec3(0, 1, 0),
					      Vec3(0, 0, 1), 
					      Vec3(-1, 0, 0)));

  RotationMatrix zyTransdMatr(combineInMatrix(Vec3(1,  0,  0),
					      Vec3(0,  0,  1), 
					      Vec3(0,  1,  0)));

  double mismatchAngle = 5.0 * M_PI / 180;

  RotationMatrix mismatchMatrix(combineInMatrix(Vec3( cos(mismatchAngle), sin(mismatchAngle), 0),
						Vec3(-sin(mismatchAngle), cos(mismatchAngle), 0), 
						Vec3(0, 0, 1)));

  RotationMatrix rm = zxTransdMatr;

  MaterialTensor tens = rotateTensor(tensMainAxis, rm);
  tens = rotateTensor(tens, mismatchMatrix);

  Vec3c unitDisplacement = Vec3c(0.0, 0.0, 1.0).norm();

  try{
    WaveMatrix waves = create_wave_matrix(n, A, f, tens, rho, unitDisplacement, FixedDisplacement);
    
    string pref = "te/x_transducer/5_deg_mismatch_";
    MatrixFFTW profile = loadFromPicture("te_cell_600x600_zy.png", 0, 1);
    MatrixFFTW profileSpectrum = profile.cleanClone();
    PlanFFTW forw(profile, profileSpectrum, FFTW_FORWARD, FFTW_ESTIMATE);
    PlanFFTW backw(profile, profileSpectrum, FFTW_BACKWARD, FFTW_ESTIMATE);

    forw.execute();
    waves.loadFftwMatrix(profileSpectrum);
    Storage spaceStor(n+1, n, n);

    for(int r=0; r<=n; r++) {
      cerr << "r = " << r << endl;

      Storage dat = waves.getStorage();
      Storage spacial = layerTransform(dat, profile, profileSpectrum, backw);
      SpacialMatrix spaceMat = getSpacialMatrixFromStorage(spacial);

      spaceMat.fillSliceWithV(r, spaceStor);

      waves.makeZShift(-dz);
    }

    MatrixFFTW sliceZX = spaceStor.sliceHW(n/2);
    MatrixFFTW sliceZY = spaceStor.sliceHD(n/2);

    saveAsPicture(sliceZX, pref + "x_zx.png");
    saveAsPicture(sliceZY, pref + "x_zy.png");
    spaceStor.save(pref+"x.txt");

  }catch(string msg){
    cout << "error: " << msg << endl;
  }
}

void volumeVisualization() {
  cerr << "volume visualization started" << endl;

  VolumeData vd(Storage(),
		Mat3(Vec3c(0, 0, -1),
		     Vec3c(0, 1, 0), 
		     Vec3c(1, 0, 0)),
		Vec3(1.3, 1.3, 1.3),
		Vec3(1.3/2, -1.3/2, -1.3/2));
	
  vd.vd.load("5_deg_mismatch_x.txt");
  cerr << endl << "volume data loaded" << endl;
  
  cerr << "0.1" << endl;
  drawVolumePovray("beam_from_x_0.1.pov", vd, 0.1);
  cerr << "0.2" << endl;
  drawVolumePovray("beam_from_x_0.2.pov", vd, 0.2);
  cerr << "0.3" << endl;
  drawVolumePovray("beam_from_x_0.3.pov", vd, 0.3);
  cerr << "0.4" << endl;
  drawVolumePovray("beam_from_x_0.4.pov", vd, 0.4);
  cerr << "0.5" << endl;
  drawVolumePovray("beam_from_x_0.5.pov", vd, 0.5);
  cerr << "0.6" << endl;
  drawVolumePovray("beam_from_x_0.6.pov", vd, 0.6);
}

void testMat3VectorProduct() {
  Mat3 m(Vec3c(0, 0, -1),
	 Vec3c(0, 1, 0), 
	 Vec3c(1, 0, 0));
    
  cout << "m = " << endl << m << endl << endl;

  cout << "m * {1, 0, 0} = " << (Vec3c(1, 0, 0) * m) << endl;
  cout << "m * {0, 1, 0} = " << (Vec3c(0, 1, 0) * m) << endl;
  cout << "m * {0, 0, 1} = " << (Vec3c(0, 0, 1) * m) << endl;

}

int main(){
  //testMat3VectorProduct();
  //workWalkOff();
  //workColorSheets();
  //dxSlowness();
  //workTellurium();
  //workTelluriumBeamStructure();
  volumeVisualization();
  //workMainTelluriumDirections();
  //wrapSplitModes();
  //testImageThroughStorage();
  //testDisplacement();
  //testNan();
  //testTransmission();
  //testHtmlBuilder();
  //testPureFourier();
  //testFreqSeq();
  //testCompositeWaveZ();
  //testRootSelectionAlgo();
  //testPresentPictures();
  //testCompositeWave();
  //testCompositeWave1();
  //testPol();
  //testPol1();
  //testStorage();
  //testStorageSaveLoad();
  //testVolumeVisualization();
  //work();
  //testSelection();
  //testStorageJson();

  // testRoots();

  // int start = (argc - 1) * 180;

  // for(int t = start + 1; t <= start + 180; ++t){
  //   cout << "start card with number " << t << endl;
  //   depthWork(t);
  // }

  //depthWork();

  //testStorageSpacialMatrix();
  //cout << "1" << endl;
}
