/**
 *  @file   LArReco/include/LArRay.h
 *
 *  @brief  Header file for LArVoxel
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_VOXEL_H
#define PANDORA_LAR_VOXEL_H 1

#include "Pandora/PandoraInputTypes.h"
#include <vector>

namespace lar_nd_reco
{

class LArVoxel
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  voxelID Total bin number for the voxel (long integer, since it can be > 2^31)
     *  @param  energyInVoxel The total deposited energy in the voxel (GeV)
     *  @param  voxelPosVect Voxel position, set as the first corner of the voxel bin
     *  @param  trackID The Geant4 ID of the (first) contributing track to this voxel
     */
    LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID);

    /**
     *  @brief  Set voxel energy (GeV)
     *
     *  @param  E voxel energy
     */
    void SetEnergy(const float E);

    /**
     *  @brief  Set voxel track id
     *
     *  @param  trackid is the track id
     */    
    void SetTrackID(const int trackid);

    long m_voxelID;                          ///< The long integer ID of the voxel (can be larger than 2^31)
    float m_energyInVoxel;                   ///< The energy in the voxel (GeV)
    pandora::CartesianVector m_voxelPosVect; ///< Position vector (x,y,z) of the first voxel corner
    int m_trackID;                           ///< The Geant4 ID of the (first) contributing track to this voxel
};

typedef std::vector<LArVoxel> LArVoxelList;

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxel::LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID) :
    m_voxelID(voxelID),
    m_energyInVoxel(energyInVoxel),
    m_voxelPosVect(voxelPosVect),
    m_trackID(trackID)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArVoxel::SetEnergy(const float E)
{
    m_energyInVoxel = E;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArVoxel::SetTrackID(const int trackid)
{
    m_trackID = trackid;
} 

//------------------------------------------------------------------------------------------------------------------------------------------

class LArVoxelProjection
{
public:
    LArVoxelProjection(const float energy, const float w, const float x, const pandora::HitType &view, const int trackid);

    float m_energy;
    float m_wire;
    float m_drift;
    pandora::HitType m_view;
    int m_trackID;
};

typedef std::vector<LArVoxelProjection> LArVoxelProjectionList;

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxelProjection::LArVoxelProjection(const float energy, const float w, const float x, 
    const pandora::HitType &view, const int trackid) :
    m_energy(energy),
    m_wire(w),
    m_drift(x),
    m_view(view),
    m_trackID(trackid)
{   
}

// namespace lar_nd_reco
}

#endif
