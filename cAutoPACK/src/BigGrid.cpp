/*
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# autopack.cpp Authors: Ludovic Autin
#
# Translation from Python initiated March 15, 2010 by Ludovic Autin
#
#
# Copyright: Graham Johnson Ludovic Autin ©2010
#
# This file "autopack.cpp" is part of autoPACK, cellPACK.
#    
#    autoPACK is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    autoPACK is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with autoPACK (See "CopyingGNUGPL" in the installation.
#    If not, see <http://www.gnu.org/licenses/>.
#
#
###############################################################################
@author: Graham Johnson, Ludovic Autin, & Michel Sanner
*/

#include <numeric>
#include <limits>
#include "BigGrid.h"
#include "BigGrid.h"
#include <algorithm>
namespace {

    openvdb::BBoxd getNeigbourBox( Ingredient * ingr, double distanceToClosestNeighbour, openvdb::Vec3d center )
    {    
        openvdb::BBoxd outerBox = ingr->getOuterBox();
        outerBox.expand(outerBox.extents()*2);
        outerBox.expand(distanceToClosestNeighbour);
        outerBox.translate(center);
        return outerBox;
    }

inline void getIJK(int u,openvdb::Coord dim,int* i_ijk){
    // = {0,0,0};    
    //openvdb::Coord ijk(0,0,0);
    int nxnynz = dim.x()*dim.y()*dim.z();
    int nynz = dim.z()*dim.y();
    //int nx = dim.x();
    int ny = dim.y();
    int nz = dim.z();
    int integer;
    double decimal;
    double fraction;
    int nysc;
    if (u < dim.z()){
        i_ijk[2] = u;
    }
    else if ((u < nynz)&&(u >= nxnynz)){
        //whats z
        fraction = (double)u/(double)dim.z();
        integer = (int) fraction;
        decimal = fraction - integer;
        i_ijk[2] = (int) round(decimal*dim.z());
        //whast y 
        i_ijk[1] = integer;  
    }
    else if ((u < nxnynz)&&(u >= nynz)){
        fraction = (double)u/(double)nynz;
        integer = (int) fraction;
        decimal = fraction - integer;
        nysc = ny * integer;
        //whast x 
        i_ijk[0] = integer;  
        fraction = (double)u/(double)nz;
        integer = (int) fraction;
        decimal = fraction - integer;
        //whats z        
        i_ijk[2] = (int) round(decimal*(double)nz);
        //whast y 
        //46867 
        //233 15477 201 77 603 77.7231
        //std::cout << integer << " " << nysc << " " << ny << " " << (int)((double)u/(double)nynz) << " " << nynz << " " << (double)u/(double)nynz << std::endl;
        i_ijk[1] = integer - (ny*(int)((double)u/(double)nynz));  
        //int (integer - (ny*int(double(u)/double(nynz))));
    }    
}

} //namespace

big_grid::big_grid( std::vector<Ingredient> const & _ingredients, double step, openvdb::Vec3d bot, openvdb::Vec3d up, unsigned seed ) :     
    distance_grid(initializeDistanceGrid(bot, up, stepsize)),
    num_points(distance_grid->activeVoxelCount()),    
    ingredientsDipatcher(_ingredients, bbox.volume(), seed),
    num_empty(bbox.volume()),
    uniform(0.0,1.0),
    half_uniform(-0.5,0.5),
    distribution(0.0,1.0),
    pickRandPt(true)
{
    std::cout << "#Grid Npoints " << distance_grid->evalActiveVoxelDim() <<  distance_grid->activeVoxelCount() << std::endl;
    generator.seed(seed);
}


openvdb::DoubleGrid::Ptr big_grid::initializeDistanceGrid( openvdb::Vec3d bot, openvdb::Vec3d up, double voxelSize )
{
    openvdb::DoubleGrid::Ptr distance_grid;
    //distance_grid = openvdb::DoubleGrid::create(dmax);
    distance_grid = openvdb::createLevelSet<openvdb::DoubleGrid>(voxelSize, spherewidth);

    //set active within the given bounding box
    // Define a coordinate with large signed indices.
    const openvdb::Vec3d ibotleft = bot;//(0,0,0)
    const openvdb::Vec3d iupright = up;//(1000,1000,10);
    openvdb::Vec3d botleft=distance_grid->worldToIndex(ibotleft);
    openvdb::Vec3d upright=distance_grid->worldToIndex(iupright);

    openvdb::Coord left(openvdb::tools::local_util::roundVec3(botleft));//(openvdb::Int32)botleft.x(),(openvdb::Int32)botleft.y(),(openvdb::Int32)botleft.z());
    openvdb::Coord right(openvdb::tools::local_util::roundVec3(upright));//(openvdb::Int32)upright.x(),(openvdb::Int32)upright.y(),(openvdb::Int32)upright.z());
    //define the active region that will be our boundary. set to max everywhere
    bbox = openvdb::CoordBBox(left,right);//min and max    
    
    const openvdb::Coord center(openvdb::tools::local_util::roundVec3(bbox.getCenter()));
    distance_grid->tree().setValueOn(center, dmax);

    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();
    std::cout << "#Testing distance access:" << std::endl;
    std::cout << "#Grid " << left << " "<< botleft << " = " << accessor_distance.getValue(left) << std::endl;
    std::cout << "#Grid " << right << " "<< upright << " = " << accessor_distance.getValue(right) << std::endl;
    std::cout << "#Grid " << bbox << std::endl;
    
    return distance_grid;
}


openvdb::Coord big_grid::getPointToDropCoord( Ingredient* ingr, double radius, double jitter, int &emptyList )
{   
    double mini_d=dmax + 1;
    openvdb::Coord mini_cijk;
    std::vector<openvdb::Coord> allIngrPts;
    
    const double dx = distance_grid->voxelSize().x();
    const double cut =  radius * dx;    
    
    if (DEBUG) std::cout << "#getPointToDropCoord " << cut << " " << mini_d <<std::endl;
    if (DEBUG) std::cout << "#retrieving available point from global grid " <<std::endl;
    for (openvdb::DoubleGrid::ValueOnIter  iter = distance_grid->beginValueOn(); iter; ++iter) {
        //for (openvdb::DoubleGrid::ValueAllIter  iter = distance_grid->beginValueAll(); iter; ++iter) {
        //before getting value check if leaf or tile    
        double d=iter.getValue();
        if (d>=cut){//the grid voxel is available and can receivethe given ingredient
            if (iter.isTileValue()){
                openvdb::CoordBBox bbox = iter.getBoundingBox();
                //iterate through it and add all the coordinates
                openvdb::Coord bbmini = bbox.min();
                openvdb::Coord bbmaxi = bbox.max();
                for (int k=bbmini.z();k<bbmaxi.z();k++){
                    for (int j=bbmini.y();j<bbmaxi.y();j++){
                        for (int i=bbmini.x();i<bbmaxi.x();i++){
                            openvdb::Coord nijk(i,j,k);
                            if (ingr->visited_rejected_coord.size() != 0 && std::find(ingr->visited_rejected_coord.begin(), ingr->visited_rejected_coord.end(), nijk) != ingr->visited_rejected_coord.end())
                            {
                                iter.setActiveState(false);
                                continue;
                            }
                            allIngrPts.push_back(nijk);
                            if (d < mini_d){                                
                                mini_d = d;
                                mini_cijk = openvdb::Coord( nijk.asVec3i() );
                            }   
                        }
                    }
                }
            } 
            else {
                openvdb::Coord cc=iter.getCoord();
                if (ingr->visited_rejected_coord.size() != 0 && (std::find(ingr->visited_rejected_coord.begin(), ingr->visited_rejected_coord.end(), cc) != ingr->visited_rejected_coord.end()))
                {
                    iter.setActiveState(false);
                    continue;
                }
                allIngrPts.push_back(cc);
                if (d < mini_d){
                    //if the point alread visisted and rejected 
                    //should be for this ingredient only
                    // not found
                    mini_d = d;
                    mini_cijk = openvdb::Coord( cc.asVec3i() );
                }
            }           
        }
    }

    
    if (DEBUG) std::cout << "#allIngrPts size " <<allIngrPts.size() << " nempty " << num_empty << " cutoff " << cut << " minid " << mini_d<<std::endl;
    if (allIngrPts.size()==0){
        std::cout << "# drop no more point \n" ;
        ingredientsDipatcher.dropIngredient(ingr); 
        emptyList = 1;
        return openvdb::Coord(0,0,0);
    }
    emptyList = 0;
    openvdb::Coord cijk;
    if (pickRandPt){
        if (ingr->packingMode=="close"){
            if (mini_d == dmax){
                cijk = getGridMiddlePoint( );
            }
            //want the smallest distance, but it is alway the same, so we get stuck here...
            //maybe use a weighting system based on distance, closed distance high prob.
            else {
                cijk = chooseTheBestPoint( allIngrPts, ingr );
            }
            
            /*Daniel - doesn't work, always came here when rejectionCounter is equal to 0; */
            if (ingr->rejectionCounter != 0 && ingr->rejectionCounter % 300 == 0){
                cijk = allIngrPts[(int)(distribution(generator) * allIngrPts.size())];
                //increase the threshold ?
                //ingr->rejectionCounter = 0;//probably not enought....risk to never end...
            }
            std::cout << "#Point to place picked:" << cijk << std::endl;
        }
        //    ptInd = 0;//sample_closest_empty(allIngrDist,allIngrPts);
        //else if (ingr->packingMode=="gradient") //&& (use_gradient)  
        //    ptInd =0;// self.gradients[ingr.gradient].pickPoint(allIngrPts) 
        //else{
        //else cijk = allIngrPts[rand() % allIngrPts.size()]; // is this working correctly?
        else {
            //or should I use a std::uniform_int_distribution<int> distribution(0,allIngrPts.size());
            cijk = allIngrPts[(int)(distribution(generator) * allIngrPts.size())];
            //rand tand to get small number first
            //cijk = allIngrPts[rand() % allIngrPts.size()];
            //why this would be on the edge and top. maybe the uniform distribution is not appropriate.
            //

        }
        //if (ptInd > allIngrPts.size()) ptInd = allIngrPts[0];            
        //}     
    }else {
        //ordered ?
        std::sort(allIngrPts.begin(),allIngrPts.end());//-(allIngrPts.size()-numActiveIngr)
        cijk = allIngrPts[0];
    }
    return cijk;
}

bool big_grid::try_drop( unsigned pid,Ingredient *ingr )
{
    //main function that decide to drop an ingredient or not
    //std::cout  <<"test_data "<< ingr.name << ' ' << pid << std::endl;        
    int i_ijk[3];
    i_ijk[0]=0;
    i_ijk[1]=0;
    i_ijk[2]=0;
    getIJK(pid, distance_grid->evalActiveVoxelDim(),i_ijk);
    openvdb::Coord cijk(i_ijk[0],i_ijk[1],i_ijk[2]);
    return try_dropCoord(cijk,ingr);
}

void big_grid::printVectorPointsToFile(const std::vector<openvdb::Coord> &allIngrPts, const std::string &fileName, const std::string &radius)
{
    std::ofstream file;
    file.open(fileName);
    for(std::size_t i=0; i<allIngrPts.size(); i++)
    {
        openvdb::Coord point = allIngrPts[i];
        file << point.x() << ' ' << point.y() << ' ' << point.z() << ' ' << radius << std::endl;
    }
    file.close();
}

openvdb::Vec3d big_grid::calculatePossition( Ingredient *ingr, openvdb::Vec3d const& offset, openvdb::math::Mat4d const& rotMatj )
{	
    openvdb::math::Transform::Ptr targetXform = openvdb::math::Transform::createLinearTransform();
    targetXform->preMult(rotMatj);
    targetXform->postTranslate(offset);
    openvdb::math::Mat4d mat = targetXform->baseMap()->getAffineMap()->getMat4();
    const openvdb::Vec3d pos = mat.transform(*(ingr->positions.begin()));
    return pos;
}

double big_grid::countDistance( std::vector<openvdb::Vec3d> const& rpossitions, Ingredient *ingr, openvdb::Vec3d const& offset, openvdb::math::Mat4d const& rotMatj )
{
    if (rpossitions.empty())
        return 0;
    const openvdb::Vec3d ingradientPosition = calculatePossition( ingr, offset, rotMatj );
    double sum = 0;
    for(auto placedIngradientsPosIter = rpossitions.cbegin(); placedIngradientsPosIter != rpossitions.cend(); placedIngradientsPosIter++)
    {
        sum += std::abs((ingradientPosition - *placedIngradientsPosIter).lengthSqr());
    }
    return sum/rpossitions.size();
}

openvdb::Vec3d color(1,0,0);
openvdb::Vec3d directions [] = {
      openvdb::Vec3d(1,0,0)
    , openvdb::Vec3d(1,1,0)
    , openvdb::Vec3d(-1,0,0)
    , openvdb::Vec3d(-1,-1,0)
    , openvdb::Vec3d(0,1,0)
    , openvdb::Vec3d(0,1,1)
    , openvdb::Vec3d(0,-1,0)
    , openvdb::Vec3d(0,-1,-1)
    , openvdb::Vec3d(0,0,1)
    , openvdb::Vec3d(0,-1,1)
    , openvdb::Vec3d(0,0,-1)
    , openvdb::Vec3d(0,1,-1)
    };

openvdb::Vec3d big_grid::findDirection(openvdb::Vec3d const& center,  double radius )
{
    const int size = sizeof(directions) / sizeof(directions[0]);
    openvdb::DoubleGrid::ConstAccessor accesor = distance_grid->getConstAccessor();
    openvdb::Vec3d retVec (1,1,1);
    double retValue = std::numeric_limits<double>::max( ) ;
    for (int i = 0; i < size; i++)
    {
        openvdb::Vec3d newPos = center + distance_grid->indexToWorld(directions[i]);
        openvdb::Coord newCoord = openvdb::Coord(openvdb::tools::local_util::roundVec3(newPos));
        const double newValue = accesor.getValue(newCoord);
        if (newValue != 0 && newValue < retValue)
        {
            retValue = newValue;
            retVec = directions[i];
        }
    }
    return retVec;
}

double big_grid::calculateValue( double i)
{   
    if (i < 0)
    {
        std::uniform_real_distribution<double> uniform (-1 ,0.0);
        double ret = uniform(generator);
        return ret;
    }
    if (i > 0)
    {
        std::uniform_real_distribution<double> uniform (0.0, 1);
        double ret = uniform(generator);
        return ret;
    }
    else
    {
        double ret =  half_uniform(generator);
        return ret;
    }
}

openvdb::Coord big_grid::findDirectionToCenter(openvdb::Coord const& point )
{    
    return getGridMiddlePoint() - point;
}

openvdb::Vec3d big_grid::generateCenterJitterOffset( openvdb::Coord const& indexCenter, openvdb::Vec3d const& ingrMaxJitter, Ingredient *ingr )
{
    openvdb::Vec3d dir = distance_grid->indexToWorld(findDirectionToCenter( indexCenter ));
    dir.normalize();
    if ( ingrMaxJitter != openvdb::Vec3d::zero() ) {
        const double x = calculateValue(dir.x());
        const double y = calculateValue(dir.y());
        const double z = calculateValue(dir.z());

        const openvdb::Vec3d randomJitter( x, y, z);
        const openvdb::Vec3d deltaOffset (ingrMaxJitter * randomJitter);

        assert( deltaOffset.lengthSqr() < ingrMaxJitter.lengthSqr() );
        return  distance_grid->indexToWorld( indexCenter + deltaOffset);
    }
    return distance_grid->indexToWorld(indexCenter) ;
}

openvdb::Vec3d big_grid::generateCloseJitterOffset( openvdb::Vec3d const& center, openvdb::Vec3d const& ingrMaxJitter, Ingredient *ingr )
{
    openvdb::Vec3d dir = findDirection(center,  1 );
    dir.normalize();
    if ( ingrMaxJitter != openvdb::Vec3d::zero() ) {
        double x = calculateValue(dir.x());
        double y = calculateValue(dir.y());
        double z = calculateValue(dir.z());

        const openvdb::Vec3d randomJitter( x, y, z);
        const openvdb::Vec3d deltaOffset (ingrMaxJitter * randomJitter);

        assert( deltaOffset.lengthSqr() < ingrMaxJitter.lengthSqr() );
        return center + distance_grid->indexToWorld( deltaOffset );
    }
    return center ;
}

openvdb::Coord big_grid::chooseTheBestPoint( const std::vector<openvdb::Coord> &allIngrPts, Ingredient *ingr )
{



    double distance = std::numeric_limits<double>::max( );
    openvdb::Coord coord (0, 0, 0);
    for ( size_t i = 0; i < allIngrPts.size(); i++ )
    {
        openvdb::Coord tempCoord = allIngrPts[i];
        const double tempDist = std::abs(countCurrentDistance( tempCoord, ingr ));
        if (tempDist < distance)
        {
            distance = tempDist;
            coord = tempCoord;
        }
    }    
    return coord;
}

openvdb::Coord big_grid::getGridMiddlePoint( )
{
    openvdb::Vec3d vec (bbox.getCenter());
    openvdb::Coord coord = openvdb::Coord(openvdb::tools::local_util::roundVec3(vec));
    return coord;
}

double big_grid::countCurrentDistance( openvdb::Coord cijk, Ingredient *ingr )
{
    openvdb::math::Mat4d rotMatj;
    rotMatj.setIdentity();
    openvdb::Vec3d center=distance_grid->indexToWorld(cijk);
    const double newDist = countDistance(rpossitions, ingr, center, rotMatj);
    return newDist;
}
 

bool big_grid::try_dropCoord( openvdb::Coord cijk,Ingredient *ingr )
{
    openvdb::Vec3d center=distance_grid->indexToWorld(cijk);

    openvdb::Vec3d globOffset;           //the point with some jitter
    openvdb::math::Mat4d globRotMatj;
    double distance = std::numeric_limits<double>::max( );
    bool placed = false;

    auto outerBox = getNeigbourBox(ingr, distance_grid->getConstAccessor().getValue(cijk), center);

    std::vector<openvdb::Vec3d> localPositions;
    std::copy_if(rpossitions.cbegin(), rpossitions.cend(), std::back_inserter(localPositions),
        [&outerBox](openvdb::Vec3d item ) { return outerBox.isInside(item);} );

    bool collision = false;
    for(unsigned i = 0; i < ingr->nbJitter; ++i) { 
        
        openvdb::Vec3d offset;
        if (collision) 
            offset = generateRandomJitterOffset(center, ingr->jitterMax);
        else
        {
            openvdb::Vec3d cc=distance_grid->worldToIndex(center);
            openvdb::Coord ci = openvdb::Coord(openvdb::tools::local_util::floorVec3(cc));
            offset = generateCenterJitterOffset(ci, ingr->jitterMax, ingr);
        }
        //const openvdb::Vec3d offset = generateRandomJitterOffset(center, ingr->jitterMax); 
        //const openvdb::Vec3d offset = generateCloseJitterOffset(center, ingr->jitterMax, ingr);
        //const openvdb::Vec3d offset = generateCenterJitterOffset(cijk, ingr->jitterMax, ingr);
        for(unsigned j = 0; j < ingr->nbJitter*3; ++j) {
            const openvdb::math::Mat4d rotMatj = generateIngredientRotation(*ingr);

            //check for collision at the given target point coordinate for the given radius     
            //collision = checkSphCollisions(offset, rotMatj, ingr->radius, ingr);
            collision = checkCollisionBasedOnGridValue(offset, rotMatj, ingr);
            if (!collision) {
                const double newDist = countDistance(localPositions, ingr, offset, rotMatj);
                if( newDist  < distance )
                {
                    if(newDist != 0)
                        distance = newDist;
                    globOffset = offset;
                    globRotMatj = rotMatj;
                    center = offset;
                    placed = true;
                    if (rtrans.empty())
                        break;
                }        
            }
        }        
    }

    if (placed)
    {
        rtrans.push_back(globOffset);
        rrot.push_back(globRotMatj);
        results.push_back(ingr);
        rpossitions.push_back( calculatePossition( ingr, globOffset, globRotMatj ) );

        placeSphereInTheGrid(globOffset, globRotMatj, ingr);
        //storePlacedIngradientInGrid(ingr, globOffset, globRotMatj);

        if (DEBUG) std::cout << "#update num_empty "<< num_empty << " " << distance_grid->activeVoxelCount() << std::endl;
    }

    return !placed;
}


bool big_grid::checkSphCollisions( openvdb::math::Vec3d const& offset,openvdb::math::Mat4d rotMatj, double radii, Ingredient* sp )
{
    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();

    //if intersect
    //should actually go through the sphere bounding box//ie object bounding 
    //value on is the shell...which could be the radius?
    //can we first test if boudning box hasOverlap ?
    for (openvdb::DoubleGrid::ValueAllCIter iter = sp->gsphere->cbeginValueAll(); iter; ++iter) {
        const openvdb::Coord sphereIndexCoord = iter.getCoord();//ijk or xyz
        const double d = iter.getValue();
        if  ( d<0 && sp->bbox.isInside(sphereIndexCoord) ) {
            //std::cout << "#sphereIndexCoord " << sphereIndexCoord << std::endl;
            openvdb::Vec3d sphereWorldCoord = sp->gsphere->indexToWorld(sphereIndexCoord);
            //std::cout << "#sphereWorldCoord " << sphereWorldCoord << " " << std::endl;
            //apply rotation
            sphereWorldCoord =rotMatj.transform(sphereWorldCoord);
            sphereWorldCoord = sphereWorldCoord + offset;
            //sphereWorldCoord = transform->indexToWorld(sphereIndexCoord);
            //sphereWorldCoord = sphereIndexCoord*matrix;
            //std::cout << "#sphereWorldCoord " << sphereWorldCoord << " " << std::endl;
            openvdb::Vec3d cc=distance_grid->worldToIndex(sphereWorldCoord);//sphereIndexCoord+woffset;//
            //std::cout << "#cc " << cc << std::endl;
            openvdb::Coord ci = openvdb::Coord(openvdb::tools::local_util::roundVec3(cc));
            //test if ci in bb ?
            //std::cout << "#ci " << ci << std::endl;
            //if  (!bbox.isInside(ci)){collide = true;return true;continue;}//or reject ?
            if  (!bbox.isInside(ci)){
                //some point are not inside the grid but could update as well?
                //continue;                       
            }//or reject ? or we can just look for collision outside ? could have an outside layer...
            //like two bounding box
            const double v  = accessor_distance.getValue(ci);
            //std::cout << "#v " << v << " " << ci << std::endl;
            //double dist = iter.getValue();
            //std::cout << "dist " << dist << std::endl;
            //check in distance if already a inside value
            if (v < 0.0) { 
                if (DEBUG) {
                    openvdb::Vec3d woffset = distance_grid->worldToIndex(offset);
                    std::cout << "#sphere position " << ci << " " << offset << " " << woffset << " reject point at d " << v <<  std::endl;
                    std::cout << "#in sphere xyz " << sp->gsphere->indexToWorld(sphereIndexCoord) << " ijk " << sphereIndexCoord << "  to grid xyz " << sphereWorldCoord << " ijk " << cc <<std::endl;
                }
                
                //reject point
                //std::cout << "reject point" << std::endl;
                //counterRej++;
                return true;
            }
        }
    }
    return false;
}

bool big_grid::checkCollisionBasedOnGridValue( openvdb::math::Vec3d const& offset,openvdb::math::Mat4d rotMatj, Ingredient* sp )
{
    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();
    for (size_t i=0; i < sp->positions.size(); i++)
    {
        openvdb::Vec3d sphereWorldCoord = sp->positions[i];
        sphereWorldCoord = rotMatj.transform(sphereWorldCoord);
        sphereWorldCoord = sphereWorldCoord + offset;

        openvdb::Vec3d cc=distance_grid->worldToIndex(sphereWorldCoord);
        openvdb::Coord ci = openvdb::Coord(openvdb::tools::local_util::floorVec3(cc));        
        const double gridValue  = accessor_distance.getValue(ci);
        const double radius = sp->radii[i];
        if ( gridValue < radius)
            return true;
    }

    return false;
}


void big_grid::placeSphereInTheGrid( openvdb::math::Vec3d const& offset,openvdb::math::Mat4d rotMatj, Ingredient* sp )
{   
    ParticeList partcieList;
    openvdb::DoubleGrid::Accessor accessor_distance = distance_grid->getAccessor();
    for (size_t i=0; i < sp->positions.size(); i++)
    {
        openvdb::Vec3R sphereWorldCoord = sp->positions[i];
        sphereWorldCoord =rotMatj.transform(sphereWorldCoord);
        sphereWorldCoord = sphereWorldCoord + offset;
        partcieList.add(sphereWorldCoord, openvdb::Real(sp->radii[i]));
    }

    openvdb::tools::ParticlesToLevelSet<openvdb::DoubleGrid, ParticeList> raster(*distance_grid);
    raster.setGrainSize(2); //a value of zero disables threading
    raster.rasterizeSpheres(partcieList);

}

openvdb::Vec3d big_grid::generateRandomJitterOffset( openvdb::Vec3d const& center, openvdb::Vec3d const& ingrMaxJitter )
{
    openvdb::Vec3d cc=distance_grid->worldToIndex(center);
    openvdb::Coord ci = openvdb::Coord(openvdb::tools::local_util::floorVec3(cc));
    openvdb::Vec3d dir = distance_grid->indexToWorld(findDirectionToCenter( ci ));
    dir.normalize();
    openvdb::Vec3d normal = dir.getArbPerpendicular();

    const auto maxJitterLength = ingrMaxJitter.lengthSqr();
    if ( maxJitterLength > 0) {
        const openvdb::Vec3d randomJitter( 
              half_uniform(generator)
            , half_uniform(generator)
            , half_uniform(generator));
        openvdb::Vec3d newIngrMaxJitter = ingrMaxJitter * 2;
        const openvdb::Vec3d deltaOffset ((newIngrMaxJitter * randomJitter) + normal);

        return center + distance_grid->indexToWorld( deltaOffset );
    }
    return center ;
}

void big_grid::storePlacedIngradientInGrid( Ingredient * ingr, openvdb::Vec3d offset, openvdb::math::Mat4d rotMatj )
{

    if (DEBUG) std::cout << "#combine the grid "<< std::endl;

    openvdb::math::Transform::Ptr targetXform =
        openvdb::math::Transform::createLinearTransform();//stepsize ?

    // Add the offset.
    const openvdb::Vec3d woffset = distance_grid->worldToIndex(offset);//offset assume the stepsize
    targetXform->preMult(rotMatj);
    targetXform->postTranslate(woffset);

    // Save copies of the two grids; compositing operations always
    // modify the A grid and leave the B grid empty.
    if (DEBUG) std::cout << "#duplicate ingredient grid "<< std::endl;

    openvdb::DoubleGrid::Ptr copyOfGridSphere = openvdb::DoubleGrid::create(dmax);            

    // Create the transformer.
    openvdb::tools::GridTransformer transformer(targetXform->baseMap()->getAffineMap()->getMat4());

    // Resample using nearest-neighbor interpolation.
    transformer.transformGrid<openvdb::tools::PointSampler, openvdb::DoubleGrid>(
        *ingr->gsphere, *copyOfGridSphere);
    copyOfGridSphere->tree().prune();
    if (DEBUG) std::cout << "#combine grid "<< std::endl;    
    openvdb::tools::csgUnion(*distance_grid, *copyOfGridSphere);

    if (DEBUG) std::cout << "#combine grid OK "<< std::endl;

    openvdb::tools::foreach(distance_grid->beginValueOn(), [this] (const openvdb::DoubleGrid::ValueOnIter& iter) 
    {
        if (!this->bbox.isInside(iter.getCoord())) iter.setActiveState(false); 
    } );
    
    num_empty = distance_grid->activeVoxelCount();    
}

openvdb::math::Mat4d big_grid::generateIngredientRotation( Ingredient const& ingr )
{
    openvdb::math::Mat4d rotMatj;
    if (ingr.useRotAxis){
        if (ingr.rotAxis.length() == 0.0)  
            rotMatj.setIdentity();
        else 
            rotMatj.setToRotation(ingr.rotAxis,uniform(generator)*ingr.rotRange);
    }else {
        rotMatj.setToRotation(openvdb::math::Vec3d(uniform(generator),uniform(generator),uniform(generator)),uniform(generator)*2.0*M_PI);
    }
    return rotMatj;
}

