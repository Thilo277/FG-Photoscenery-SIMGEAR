///@file
/// Sound effect management class
///
/// Provides a sound manager class to keep track of multiple sounds and manage
/// playing them with different effects and timings.
//
// Sound manager initially written by David Findlay
// <david_j_findlay@yahoo.com.au> 2001
//
// C++-ified by Curtis Olson, started March 2001.
// Modified for the new SoundSystem by Erik Hofman, October 2009
//
// Copyright (C) 2001  Curtis L. Olson - http://www.flightgear.org/~curt
// Copyright (C) 2009-2019 Erik Hofman <erik@ehofman.com>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _SG_SOUNDMGR_HXX
#define _SG_SOUNDMGR_HXX 1

#include <string>
#include <vector>
#include <map>
#include <memory> // for std::unique_ptr

#include <simgear/compiler.h>
#include <simgear/structure/subsystem_mgr.hxx>
#include <simgear/math/SGMath.hxx>


// Speed of sound in meters per second
#define SPEED_OF_SOUND		340.3

// forward decls
class SGSampleGroup;
class SGSoundSample;

/**
 * Manage a collection of SGSampleGroup instances
 */
class SGSoundMgr : public SGSubsystem
{
public:
    SGSoundMgr();
    ~SGSoundMgr();

    // Subsystem API.
    void init() override;
    void reinit() override;
    void resume() override;
    void suspend() override;
    void update(double dt) override;

    // Subsystem identification.
    static const char* staticSubsystemClassId() { return "sound"; }

    void stop();

    /**
     * Select a specific sound device.
     * Requires a init/reinit call before sound is actually switched.
     */
    inline void select_device(const char* devname) { _device_name = devname; }

    /**
     * Test is the sound manager is in a working condition.
     * @return true is the sound manager is working
     */
    bool is_working() const;

    /**
     * Set the sound manager to a  working condition.
     */
    void activate();

    /**
     * Test is the sound manager is in an active and working condition.
     * @return true is the sound manager is active
     */
    inline bool is_active() const { return _active; }

    /**
     * Register a sample group to the sound manager.
     * @param sgrp Pointer to a sample group to add
     * @param refname Reference name of the sample group
     * @return true if successful, false otherwise
     */
    bool add( SGSampleGroup *sgrp, const std::string& refname );

    /**
     * Remove a sample group from the sound manager.
     * @param refname Reference name of the sample group to remove
     * @return true if successful, false otherwise
     */
    bool remove( const std::string& refname );

    /**
     * Test if a specified sample group is registered at the sound manager
     * @param refname Reference name of the sample group test for
     * @return true if the specified sample group exists
     */
    bool exists( const std::string& refname );

    /**
     * Find a specified sample group in the sound manager
     *
     * @param refname Reference name of the sample group to find
     * @param create  If the group should be create if it does not exist
     * @return A pointer to the SGSampleGroup
     */
    SGSampleGroup *find( const std::string& refname, bool create = false );

    /**
     * Set the Cartesian position of the sound manager.
     *
     * @param pos OpenAL listener position
     */
    void set_position( const SGVec3d& pos, const SGGeod& pos_geod );

    /**
     * Get the position of the sound manager.
     * This is in the same coordinate system as OpenGL; y=up, z=back, x=right
     *
     * @return OpenAL listener position
     */
    const SGVec3d& get_position() const;

    /**
     * Set the velocity vector (in meters per second) of the sound manager
     * This is the horizontal local frame; x=north, y=east, z=down
     *
     * @param vel Velocity vector
     */
    void set_velocity( const SGVec3d& vel ) {
        _velocity = vel; _changed = true;
    }

    /**
     * Get the velocity vector of the sound manager
     * This is in the same coordinate system as OpenGL; y=up, z=back, x=right.
     *
     * @return Velocity vector of the OpenAL listener
     */
    inline SGVec3f get_velocity() { return toVec3f(_velocity); }

    /**
     * Set the orientation of the sound manager
     *
     * @param ori Quaternation containing the orientation information
     */
    void set_orientation( const SGQuatd& ori );

    /**
     * Get the orientation of the sound manager
     *
     * @return Quaternation containing the orientation information
     */
    const SGQuatd& get_orientation() const;

    /**
     * Get the direction vector of the sound manager
     * This is in the same coordinate system as OpenGL; y=up, z=back, x=right.
     *
     * @return Look-at direction of the OpenAL listener
     */
    SGVec3f get_direction() const;

    enum {
        NO_SOURCE = (unsigned int)-1,
        NO_BUFFER = (unsigned int)-1,
        FAILED_BUFFER = (unsigned int)-2
    };

    /**
     * Set the master volume.
     *
     * @param vol Volume (must be between 0.0 and 1.0)
     */
    void set_volume( float vol );

    /**
     * Get the master volume.
     *
     * @return Volume (must be between 0.0 and 1.0)
     */
    inline float get_volume() { return _volume; }

    /**
     * Set the speed of sound.
     *
     * @param vel Sound velocity
     */
    void set_sound_velocity( double vel ) { _sound_velocity = vel; }

    /**
     * Get a free OpenAL source-id
     *
     * @return NO_SOURCE if no source is available
     */
    unsigned int request_source();

    /**
     * Free an OpenAL source-id for future use
     *
     * @param source OpenAL source-id to free
     */
    void release_source( unsigned int source );

    /**
     * Get a free OpenAL buffer-id
     * The buffer-id will be assigned to the sample by calling this function.
     *
     * @param sample Pointer to an audio sample to assign the buffer-id to
     * @return NO_BUFFER if loading of the buffer failed.
     */
    unsigned int request_buffer(SGSoundSample *sample);

    /**
     * Free an OpenAL buffer-id for this sample
     *
     * @param sample Pointer to an audio sample for which to free the buffer
     */
    void release_buffer( SGSoundSample *sample );

    /**
     * Initialize sample for playback.
     *
     * @param sample Pointer to an audio sample to initialize.
     */
    void sample_init( SGSoundSample *sample );

    /**
     * Stop and destroy a sample
     *
     * @param sample Pointer to an audio sample to destroy.
     */
    void sample_destroy( SGSoundSample *sample );

    /**
     * Start playback of a sample
     *
     * @param sample Pointer to an audio sample to start playing.
     */
    void sample_play( SGSoundSample *sample );

    /**
     * Stop a sample
     *
     * @param sample Pointer to an audio sample to stop.
     */
    void sample_stop( SGSoundSample *sample );

    /**
     * Suspend playback of a sample
     *
     * @param sample Pointer to an audio sample to suspend.
     */
    void sample_suspend( SGSoundSample *sample );

    /**
     * Resume playback of a sample
     *
     * @param sample Pointer to an audio sample to resume.
     */
    void sample_resume( SGSoundSample *sample );

    /**
     * Check if a sample is stopped, or still playing
     *
     * @param sample Pointer to an audio sample to test.
     * @return true if the sample is stopped.
     */
    bool is_sample_stopped( SGSoundSample *sample );

    /**
     * Update all status and 3d parameters of a sample.
     *
     * @param sample Pointer to an audio sample to update.
     */
    void update_sample_config( SGSoundSample *sample, SGVec3d& position, SGVec3f& orientation, SGVec3f& velocity );

    /**
     * Test if the position of the sound manager has changed.
     * The value will be set to false upon the next call to update_late()
     *
     * @return true if the position has changed
     */
    inline bool has_changed() { return _changed; }

    /**
     * Load a sample file and return it's configuration and data.
     *
     * @param samplepath Path to the file to load
     * @param data Pointer to a variable that points to the allocated data
     * @param format Pointer to a vairable that gets the OpenAL format
     * @param size Pointer to a vairable that gets the sample size in bytes
     * @param freq Pointer to a vairable that gets the sample frequency in Herz
     * @return true if succesful, false on error
     */
    virtual bool load( const std::string &samplepath,
                       void** data,
                       int *format,
                       size_t *size,
                       int *freq,
                       int *block );

    /**
     * Get a list of available playback devices.
     */
    std::vector<std::string> get_available_devices();

    /**
     * Get the current OpenAL vendor or rendering backend.
     */
    const std::string& get_vendor() { return _vendor; }
    const std::string& get_renderer() { return _renderer; }

    bool testForError(std::string s, std::string name = "sound manager");

private:
    class SoundManagerPrivate;
    /// private implementation object
    std::unique_ptr<SoundManagerPrivate> d;

    bool _block_support;
    bool _active = false;
    bool _changed = true;
    float _volume = 0.0f;

    double _sound_velocity = SPEED_OF_SOUND;

    // Position of the listener.
    SGGeod _geod_pos;

    // Velocity of the listener.
    SGVec3d _velocity = SGVec3d::zeros();

    std::string _renderer = "unknown";
    std::string _vendor = "unknown";
    std::string _device_name;

    bool testForALCError(std::string s);
    bool testForError(void *p, std::string s);

    void update_sample_config( SGSampleGroup *sound );
};


#endif // _SG_SOUNDMGR_HXX
