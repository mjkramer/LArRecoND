/**
 *  @file   LArReco/include/LArNDLArGeomSimple.h
 *
 *  @brief  Header file for storing geometry information for ND LAr
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_ND_LAR_GEOM_SIMPLE_H
#define PANDORA_LAR_ND_LAR_GEOM_SIMPLE_H 1

#include "Pandora/PandoraInputTypes.h"

namespace lar_nd_reco
{

class LArNDLArTPCSimple
{
public:

    LArNDLArTPCSimple();

    LArNDLArTPCSimple(const float x_min, const float x_max, const float y_min, const float y_max,
                 const float z_min, const float z_max, const unsigned int tpcID);

    bool IsInTPC(const pandora::CartesianVector &pos) const;

    float m_x_min;
    float m_x_max;
    float m_y_min;
    float m_y_max;
    float m_z_min;
    float m_z_max;
    unsigned int m_TPC_ID;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNDLArTPCSimple::LArNDLArTPCSimple() :
    m_x_min{0.f}, m_x_max{0.f},
    m_y_min{0.f}, m_y_max{0.f},
    m_z_min{0.f}, m_z_max{0.f},
    m_TPC_ID{999}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNDLArTPCSimple::LArNDLArTPCSimple(const float x_min, const float x_max, const float y_min, const float y_max,
                                  const float z_min, const float z_max, const unsigned int tpcID) :
    m_x_min{x_min}, m_x_max{x_max},
    m_y_min{y_min}, m_y_max{y_max},
    m_z_min{z_min}, m_z_max{z_max},
    m_TPC_ID{tpcID}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArNDLArTPCSimple::IsInTPC(const pandora::CartesianVector &pos) const
{
    if (pos.GetX() < m_x_min || pos.GetX() > m_x_max)
        return false;
    if (pos.GetY() < m_y_min || pos.GetY() > m_y_max)
        return false;
    if (pos.GetZ() < m_z_min || pos.GetZ() > m_z_max)
        return false;

    return true;
}

class LArNDLArGeomSimple
{
public:

    /**
     *  @brief  default constructor
     */  
    LArNDLArGeomSimple();

    /**
     *  @brief  Get the TPC number from a 3D position
     *
     *  @param  position 3d position to query
     *
     *  @return tpc number as unsigned int
     */
    unsigned int GetTPCNumber(const pandora::CartesianVector &position) const;

    /**
     *  @brief  Get the module number from a 3D position
     *
     *  @param  position 3d position to query
     *
     *  @return module number as unsigned int
     */
    unsigned int GetModuleNumber(const pandora::CartesianVector &position) const;

    // I'm going to label each TPC from 0 -> 13 for the first row (low z)
    // and then 15 -> 26 for the second row etc
    std::map<unsigned int,LArNDLArTPCSimple> m_TPCs;

private:

    void PopulateBoundaries();
};


//------------------------------------------------------------------------------------------------------------------------------------------

LArNDLArGeomSimple::LArNDLArGeomSimple()
{
    this->PopulateBoundaries();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArNDLArGeomSimple::GetTPCNumber(const pandora::CartesianVector &position) const
{
    for (auto const &tpc : m_TPCs)
    {
        if (tpc.second.IsInTPC(position))
            return tpc.first;
    }
    return 999;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArNDLArGeomSimple::GetModuleNumber(const pandora::CartesianVector &position) const
{
    for (auto const &tpc : m_TPCs)
    {
        if (tpc.second.IsInTPC(position))
            return tpc.first / 2;
    }
    return 999;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNDLArGeomSimple::PopulateBoundaries()
{
    // These numbers are roughly accurate to 0.1mm and obtained from looking at hit positions
    const float m_x_min = -347.85; 
//    const float m_x_max = 349.7;
    const float m_y_min = -216.7; 
    const float m_y_max = 83.0;
    const float m_z_min = 417.92;
//    const float m_z_max = 913.6; 

    // The repeating pattern size (same in x and z)
    const float repeat_size = 100.0;

    // Gap sizes
    const float x_cathode_gap = 0.62;
    const float x_tpc_width = 47.54;
    const float z_module_gap = 4.33;

    // These are hardcoded numbers for now as we don't have anything better available
    const unsigned int m_n_modules_x = 7;
    const unsigned int m_n_modules_z = 5;

    for (unsigned int mx = 0; mx < m_n_modules_x; ++mx)
    {
        for (unsigned int mz = 0; mz < m_n_modules_z; ++mz)
        {
            // Each module contains two TPCs: left (low x) and right (high x)
            float x_min{m_x_min + mx*repeat_size};
            float x_max{x_min + x_tpc_width};
            const float z_min{m_z_min + mz*repeat_size};
            const float z_max{z_min + repeat_size - z_module_gap};
            const unsigned int moduleID = mx + (mz * m_n_modules_x);
            LArNDLArTPCSimple leftTPC(x_min,x_max,m_y_min,m_y_max,z_min,z_max,2*moduleID);

            // Now for the right TPC
            x_min = x_max + x_cathode_gap;
            x_max = x_min + x_tpc_width;
            LArNDLArTPCSimple rightTPC(x_min,x_max,m_y_min,m_y_max,z_min,z_max,2*moduleID+1);

            // Store the TPCs in the map, keyed on their tpcID
            m_TPCs[leftTPC.m_TPC_ID] = leftTPC;
            m_TPCs[rightTPC.m_TPC_ID] = rightTPC;            
        }
    }
}

} // namespace lar_nd_reco

#endif
