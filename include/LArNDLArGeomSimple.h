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

    /**
     *  @brief  default constructor
     */  
    LArNDLArTPCSimple();

    /**
     *  @brief  constructor with coordinate limits and id
     *
     *  @param  x_min minimum x value of the TPC box
     *  @param  x_max maximum x value of the TPC box
     *  @param  y_min minimum y value of the TPC box
     *  @param  y_max maximum y value of the TPC box
     *  @param  z_min minimum z value of the TPC box
     *  @param  z_max maximum z value of the TPC box
     *  @param  tpcID unique id of the tpc
     */  
    LArNDLArTPCSimple(const double x_min, const double x_max, const double y_min, const double y_max,
                 const double z_min, const double z_max, const unsigned int tpcID);

    bool IsInTPC(const pandora::CartesianVector &pos) const;

    double m_x_min;         ///< minimum x value of the TPC cubioid
    double m_x_max;         ///< maximum x value of the TPC cubioid
    double m_y_min;         ///< minimum y value of the TPC cubioid
    double m_y_max;         ///< maximum y value of the TPC cubioid
    double m_z_min;         ///< minimum z value of the TPC cubioid
    double m_z_max;         ///< maximum z value of the TPC cubioid
    unsigned int m_TPC_ID;  ///< unique id of the TPC

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

inline LArNDLArTPCSimple::LArNDLArTPCSimple(const double x_min, const double x_max, const double y_min, const double y_max,
                                  const double z_min, const double z_max, const unsigned int tpcID) :
    m_x_min{x_min}, m_x_max{x_max},
    m_y_min{y_min}, m_y_max{y_max},
    m_z_min{z_min}, m_z_max{z_max},
    m_TPC_ID{tpcID}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArNDLArTPCSimple::IsInTPC(const pandora::CartesianVector &pos) const
{
    if (pos.GetX() <= m_x_min || pos.GetX() >= m_x_max)
        return false;
    if (pos.GetY() <= m_y_min || pos.GetY() >= m_y_max)
        return false;
    if (pos.GetZ() <= m_z_min || pos.GetZ() >= m_z_max)
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

    /**
     *  @brief  Add a TPC to the geometry
     *
     *  @param  min_x minimum x position of the TPC
     *  @param  max_x maximum x position of the TPC
     *  @param  min_y minimum y position of the TPC
     *  @param  max_y maximum y position of the TPC
     *  @param  min_z minimum z position of the TPC
     *  @param  max_z maximum z position of the TPC
     *  @param  tpcID tpc id as unsigned int
     */
    void AddTPC(const double min_x, const double max_x, const double min_y, const double max_y,
                const double min_z, const double max_z, const unsigned int tpcID);

    /**
     *  @brief  Get the box surrounding all TPCs
     *
     *  @param  min_x minimum x position of the TPCs
     *  @param  max_x maximum x position of the TPCs
     *  @param  min_y minimum y position of the TPCs
     *  @param  max_y maximum y position of the TPCs
     *  @param  min_z minimum z position of the TPCs
     *  @param  max_z maximum z position of the TPCs
     */
    void GetSurroundingBox(double &min_x, double &max_x, double &min_y, double &max_y, double &min_z, double &max_z) const;

    std::map<unsigned int,LArNDLArTPCSimple> m_TPCs; ///< map of the TPCs keyed on their unique id

private:

};


//------------------------------------------------------------------------------------------------------------------------------------------

LArNDLArGeomSimple::LArNDLArGeomSimple()
{
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

inline void LArNDLArGeomSimple::AddTPC(const double min_x, const double max_x, const double min_y, const double max_y,
                   const double min_z, const double max_z, const unsigned int tpcID)
{
    if (m_TPCs.count(tpcID))
        std::cout << "LArNDLArGeomSimple: trying to add another TPC with tpc id " << tpcID << "! Doing nothing. " << std::endl;
    else
        m_TPCs[tpcID] = LArNDLArTPCSimple(min_x, max_x ,min_y, max_y, min_z, max_z, tpcID);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNDLArGeomSimple::GetSurroundingBox(double &min_x, double &max_x, double &min_y, double &max_y, double &min_z, double &max_z) const
{
    min_x = std::numeric_limits<double>::max();
    min_y = std::numeric_limits<double>::max();
    min_z = std::numeric_limits<double>::max();
    max_x = std::numeric_limits<double>::min();
    max_y = std::numeric_limits<double>::min();
    max_z = std::numeric_limits<double>::min();

    for (auto const &tpc : m_TPCs)
    {
        min_x = tpc.second.m_x_min < min_x ? tpc.second.m_x_min : min_x; 
        min_y = tpc.second.m_y_min < min_y ? tpc.second.m_y_min : min_y; 
        min_z = tpc.second.m_z_min < min_z ? tpc.second.m_z_min : min_z; 
        max_x = tpc.second.m_x_max > max_x ? tpc.second.m_x_max : max_x; 
        max_y = tpc.second.m_y_max > max_y ? tpc.second.m_y_max : max_y; 
        max_z = tpc.second.m_z_max > max_z ? tpc.second.m_z_max : max_z; 
    }
}

} // namespace lar_nd_reco

#endif
