/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    faceAgglomerate

Description

    Agglomerate boundary faces using the pairPatchAgglomeration algorithm.
    It writes a map from the fine to coarse grid.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "CompactListList.H"
#include "unitConversion.H"
#include "labelListIOList.H"
#include "syncTools.H"
#include "globalIndex.H"
#include <assert.h>


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
// In openFoam for each region one viewFactor Matrix will be created:
// This function is consisted of four parts: 
// a. For all the specified regions, check the viewFactorDict and calculate three variables: 
// 	1. Total faces in each region
// 	2. Total faces in each region that the user has set the "nFacesInCoarsestLevel" parameter for them in viewFactorsDict
//	3. An approximation on the number of the coarsed faces in each region for faces that the user has set the
//		"nFacesInCoarsestLevel" parameter for them in viewFactorsDict
//
// b. Total amount of available memory for viewFactor matrixes of all the regions
//	will be acquired from the amount that the user has specified in controlDict
//
// c. some evaluation will be performed to calculate the "nFacesInCoarsestLevel" parameter for the faces that the user 
//	has not specified.
//
// d. write the calculated value of nFacesInCoarsestLevel for rest of faces in the viewFactorDict of each region 


  int argc_aux = 3;
  char** argv_aux;
  argv_aux = new char*[3];
  argv_aux[0]="prefaceAgglomerate";
  argv_aux[1]="-region";
  argv_aux[2]=argv[1];


  #include "addRegionOption.H"
  #include "addDictOption.H"
	
  //
  // setRootCase.H
  // ~~~~~~~~~~~~~
  Foam::argList args(argc_aux, argv_aux);
  if (!args.checkRootCase())
  {
	Foam::FatalError.exit();
  }

  #include "createTime.H"


  const int numRegions = argc-1;
  List<word> regions(numRegions);
  for(int i=0; i<numRegions; i++) regions[i]=argv[i+1];
  List<scalar> region_totalFaces(numRegions);
  List<scalar> region_totalFacs_With_viewFactorsDict_specified(numRegions);
  List<scalar> region_Approx_nCoarsed_Faces(numRegions);

//-----------------------------------------------------------------------------------------------------
// a. This part performs the estimations on the faces in each region
//-----------------------------------------------------------------------------------------------------

  forAll(regions, id)
  {
      	Info << "region: " << regions[id] << endl;

	Foam::word regionName = regions[id];
	// create Named Mesh for this region
	// ~~~~~~~~~~~~~~~~~
	    Foam::fvMesh mesh
	    (
		Foam::IOobject
		(
		    regionName,
		    runTime.timeName(),
		    runTime,
		    Foam::IOobject::MUST_READ
		)
	    );


   	 word dictName("viewFactorsDict");
	// read the view factor dictionary for this region
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    IOobject dictIO
	    (
		dictName,
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	    );

	IOdictionary agglomDict(dictIO);

	const polyBoundaryMesh& boundary = mesh.boundaryMesh();

	label totalFaces = 0;
	label totalFacs_With_viewFactorsDict_specified = 0;
	label Approx_nCoarsed_Faces = 0;


	forAll(boundary, patchId) { totalFaces += boundary[patchId].size(); }

	forAllConstIter(dictionary, agglomDict, iter)
	{
	    labelList patchIds = boundary.findIndices(iter().keyword());


            forAll(patchIds, i)
            {
            	label patchI =  patchIds[i];
	        const polyPatch& pp = boundary[patchI];
   		Info << "patch name: " << pp.name() << endl;
		Info << "patch size: " << pp.size() << endl;

 	    	totalFacs_With_viewFactorsDict_specified += boundary[patchI].size();

 		if (!pp.coupled())
	        {
	           label nFacesInCoarsestLevel_ = readLabel(agglomDict.subDict(pp.name()).lookup("nFacesInCoarsestLevel"));
    	           Approx_nCoarsed_Faces += boundary[patchI].size()/nFacesInCoarsestLevel_;
                }
            }
	}

    	Info << "Total number of faces in this region is: " << totalFaces << endl;
	Info << "Total number of facs with viewFactorsDict specified by user: " << totalFacs_With_viewFactorsDict_specified << endl;
    	Info << "Approxate number of Coarsed Faces based on the parameters that the user specified: " << Approx_nCoarsed_Faces << endl;

	region_totalFaces[id] = totalFaces;
	region_totalFacs_With_viewFactorsDict_specified[id] = totalFacs_With_viewFactorsDict_specified;
	region_Approx_nCoarsed_Faces[id] = Approx_nCoarsed_Faces;

  }

//-----------------------------------------------------------------------------------------------------
// b. This part reads the total amount of available memory for viewFactor matrixes of all the regions
//-----------------------------------------------------------------------------------------------------

  // read the available memory for viewfactor calculations
  const word dictName("controlDict");
  IOdictionary controlDict
  (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
   );

   if (!controlDict.found("MaxMemoryForViewFactors")) {
        Info << "The maximum memory allocated for viewFactor is not specified. \nUsing the defaul value: 1000MB " << endl;
   }
   const double totalMem = controlDict.lookupOrDefault<scalar>("MaxMemoryForViewFactors", 1000) * 1000000;
   Info << "Total allocated memory for view factor matrices is: " << totalMem << " bytes" << endl;


//-----------------------------------------------------------------------------------------------------
// c. This part performs some calculations to evaluate the "nFacesInCoarsestLevel" parameter for the faces that the user 
//	has not specified
//-----------------------------------------------------------------------------------------------------
   double occupiedMemory_by_faces_With_viewFactorsDict_specified = 0.0;
   forAll(regions, id)
   {
     	occupiedMemory_by_faces_With_viewFactorsDict_specified += Foam::pow( region_Approx_nCoarsed_Faces[id],2.0 ) * sizeof(double);
   }
   assert( (occupiedMemory_by_faces_With_viewFactorsDict_specified<totalMem)
		 && "The faces that the user specified thir parameters occupy more memory than the amount allocated in controlDict!\n");


  label nFacesInCoarsestLevel_rest_of_faces = 0;
  double occupiedMemory = totalMem+1.0;
  while( occupiedMemory>totalMem )
  {
	++nFacesInCoarsestLevel_rest_of_faces;
	occupiedMemory = 0.0;
	forAll(regions, id)
	{
	    double nfaces = ( region_totalFaces[id]-region_totalFacs_With_viewFactorsDict_specified[id] ) /
						 nFacesInCoarsestLevel_rest_of_faces + region_Approx_nCoarsed_Faces[id];
	    occupiedMemory += Foam::pow( nfaces,2.0 ) * sizeof(double);
	}
  }

//-----------------------------------------------------------------------------------------------------
// d. write the calculated value of nFacesInCoarsestLevel for rest of faces in the viewFactorDict of each region 
//-----------------------------------------------------------------------------------------------------

  forAll(regions, id)
  {
	Foam::word regionName = regions[id];
	// create Named Mesh for this region
	// ~~~~~~~~~~~~~~~~~
	    Foam::fvMesh mesh
	    (
		Foam::IOobject
		(
		    regionName,
		    runTime.timeName(),
		    runTime,
		    Foam::IOobject::MUST_READ
		)
	    );

   	word dictName("viewFactorsDict");
	IOobject dictIO
	(
	    dictName,
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ_IF_MODIFIED,
	    IOobject::NO_WRITE
	);

	// Read control dictionary
	IOdictionary agglomDict(dictIO);

	dictionary& rest_of_faces_Dict = agglomDict.subDict("rest_of_faces");
	rest_of_faces_Dict.set("nFacesInCoarsestLevel",nFacesInCoarsestLevel_rest_of_faces);
	rest_of_faces_Dict.set("featureAngle",10);

	agglomDict.regIOobject::write();
   }

   Info<< "End\n" << endl;
   return 0;
}


// ************************************************************************* //
