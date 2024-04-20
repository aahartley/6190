//-------------------------------------------------------
//
//  DynamicalState.h
//
//  Container for data associated with the dynamics
//  degrees of freedom in a system.
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//--------------------------------------------------------

#ifndef ____PBA_DYNAMICALSTATE_H____
#define ____PBA_DYNAMICALSTATE_H____

#include "Vector.h"
#include "Color.h"
#include "AABB.h"
#include <vector>
#include <string>
#include <map>
#include <memory>


namespace pba
{

//! This class handles data for a single attribute
template<typename T>
class DSAttribute
{
  public:

    DSAttribute() : name("unknown") {}
    DSAttribute( const std::string& nam, const T& def ) : name(nam), defVal(def) {}
   ~DSAttribute(){}

    const size_t size() const { return data.size(); }
    const bool empty() const { return data.empty(); }
    void set(size_t i, const T& value ) { data[i] = value; }
    const T& get(size_t i ) const { return data[i]; }
    T& get(size_t i ) { return data[i]; }
    void expand_to( size_t n )
    {
       if( data.size() >= n ){ return; }
       size_t old_size = data.size();
       data.resize(n);
       for( size_t i=old_size;i<data.size();i++ )
       {
          data[i] = defVal;
       }
    }
    void clear() { data.clear(); }
    const std::string& attr_name() const { return name; }
    const T& default_value() const { return defVal; }
    typename std::vector<T>::const_iterator cbegin() const { return data.begin(); }
    typename std::vector<T>::const_iterator cend() const { return data.end(); }
    typename std::vector<T>::iterator begin() { return data.begin(); }
    typename std::vector<T>::iterator end() { return data.end(); }

    void erase(const size_t i)
    {
        typename std::vector<T>::iterator removable = data.begin() + i;
	data.erase(removable);
    }

  private:
    std::vector<T> data;
    std::string name;
    T defVal;
};

//! This class is a collection of attributes for many particles. Some default attributes exist, and more can be created.
/*!
     Any number of attributes can be created.  Attributes must have int, float, Vector, or Color type.  
     There are a set of pre-built attributes.  They are
     Name                  Type
     pos (position)        Vector
     vel (velocity)        Vector
     accel (acceleration)  Vector
     mass                  float
     id                    int
     ci                    Color
 */
class DynamicalStateData 
{
  public:

    DynamicalStateData( const std::string& nam = "DynamicDataNoName");
    DynamicalStateData( const DynamicalStateData& d );
   ~DynamicalStateData();

    DynamicalStateData& operator= ( const DynamicalStateData& d );

    //! Create a new attribute of type int
    void create_attr( const std::string& nam, const int& def );
    //! Create a new attribute of type float
    void create_attr( const std::string& nam, const float& def );
    //! Create a new attribute of type Vector
    void create_attr( const std::string& nam, const Vector& def );
    //! Create a new attribute of type Color
    void create_attr( const std::string& nam, const Color& def );

    // Add a single particle
    const size_t add();
    // Add many particles
    const size_t add( const size_t nb );
    //! Return the number of particles
    size_t nb() const;
    //! Remove all particles
    void clear();

    //! Return the int attribute named, for the particle indicated
    const int& get_int_attr( const std::string& nam, const size_t p ) const;
    //! Return the float attribute named, for the particle indicated
    const float& get_float_attr( const std::string& nam, const size_t p ) const;
    //! Return the Vector attribute named, for the particle indicated
    const Vector& get_vector_attr( const std::string& nam, const size_t p ) const;
    //! Return the Color attribute named, for the particle indicated
    const Color& get_color_attr( const std::string& nam, const size_t p ) const;
   
    //! Return the pos (position) value for the particle indicated 
    const Vector& pos( const size_t p ) const;
    //! Return the velocity value for the particle indicated 
    const Vector& vel( const size_t p ) const;
    //! Return the acceleration value for the particle indicated 
    const Vector& accel( const size_t p ) const;
    //! Return the mass value for the particle indicated 
    const float& mass( const size_t p ) const;
      //! Return the mass value for the particle indicated 
    const float& rad( const size_t p ) const;
    //! Return the id value for the particle indicated 
    const int& id( const size_t p ) const;
    //! Return the Color value for the particle indicated 
    const Color& ci( const size_t p ) const;

    //! Set the value of the named int attribute and designated particle 
    void set_attr( const std::string& nam, const size_t p, const int& value ); 
    //! Set the value of the named float attribute and designated particle 
    void set_attr( const std::string& nam, const size_t p, const float& value ); 
    //! Set the value of the named Vector attribute and designated particle 
    void set_attr( const std::string& nam, const size_t p, const Vector& value ); 
    //! Set the value of the named Color attribute and designated particle 
    void set_attr( const std::string& nam, const size_t p, const Color& value ); 

    //! Set the pos attribute for the designated particle
    void set_pos( const size_t p, const Vector& value );
    //! Set the vel attribute for the designated particle
    void set_vel( const size_t p, const Vector& value );
    //! Set the accel attribute for the designated particle
    void set_accel( const size_t p, const Vector& value );
    //! Set the mass attribute for the designated particle
    void set_mass( const size_t p, const float& value );
      //! Set the mass attribute for the designated particle
    void set_particle_radius( const size_t p, const float& value );
    //! Set the id attribute for the designated particle
    void set_id( const size_t p, const int& value );
    //! Set the ci attribute for the designated particle
    void set_ci( const size_t p, const Color& value );

    //! Returns the number of particles ersased
    int erase_outside_bounds( const Vector& llc, const Vector& urc );

    //! Return a list of names of int attributes
    std::vector<std::string> show_int_attrs() const;
    //! Return a list of names of float attributes
    std::vector<std::string> show_float_attrs() const;
    //! Return a list of names of Vector attributes
    std::vector<std::string> show_vector_attrs() const;
    //! Return a list of names of Color attributes
    std::vector<std::string> show_color_attrs() const;
    //! Return a list of all of the attributes
    std::vector<std::string> show_all_attrs() const;

    //! Return whether an attribute exists with the designated name
    bool attr_exists( const std::string& nam ) const;

    //! Merge the designated state data with this one
    void merge( const DynamicalStateData& g );

    //! Name assigned to this object
    const std::string& Name() const { return name; }

    //! Current absolute simulation tiem
    const double time() const { return t; }
    //! Increment the absolute time
    void update_time(const double dt){ t += dt; }

  protected:

    double t;
    size_t nb_items;

    // attribute collections
    std::map< std::string, DSAttribute<int>  > int_attributes;
    std::map< std::string, DSAttribute<float>  > float_attributes;
    std::map< std::string, DSAttribute<Vector> > vector_attributes;  
    std::map< std::string, DSAttribute<Color> > color_attributes;  

    std::map< std::string, DSAttribute<Vector> >::iterator    positions;
    std::map< std::string, DSAttribute<Vector> >::iterator    velocities;
    std::map< std::string, DSAttribute<Vector> >::iterator    accelerations;
    std::map< std::string, DSAttribute<float> >::iterator     masses;
    std::map< std::string, DSAttribute<float> >::iterator     radiuses;
    std::map< std::string, DSAttribute<int> >::iterator       ids;
    std::map< std::string, DSAttribute<Color> >::iterator     cis;

    std::string name;

    void re_find_main_attrs();


};


//! DynamicalState is a smart pointer of DynamicalStateData
typedef std::shared_ptr<DynamicalStateData> DynamicalState;

//! Create a smart pointer for a dynamical state
DynamicalState CreateDynamicalState( const std::string& nam = "DynamicalDataNoName" );

//! Make a deep copy of a DynamicalState
DynamicalState copy( const DynamicalState d );

//! Average of the pos attribute
Vector geometric_center( const DynamicalState& d );
//! Bounding box enclosing all pos attribute locations
AABB BoundingBox( const DynamicalState& d );
AABB BoundingBox( const DynamicalStateData& d );
//! Global translate of all positions 
void Translate( const Vector& delta, DynamicalState& d );
//! Global scaling of all positions 
void Scale( const Vector& reference, const Vector& scale, DynamicalState& d );
//! Global rotation of all positions 
void Rotate( const Vector& reference, const Vector& rotor, DynamicalState& d );


}
#endif
