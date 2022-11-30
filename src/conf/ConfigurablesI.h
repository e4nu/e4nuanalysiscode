/**
 * \info These parameters are configurable 
 * the default values are set here
 **/

#ifndef _CONFIGURABLES_I_H_
#define _CONFIGURABLES_I_H_

namespace e4nu { 

  namespace conf {

    class e4nuConfigurationI {

    public: 
      // Default constructor
      e4nuConfigurationI() ;

      // Configure with input file 
      e4nuConfigurationI( std::string input_file ) ;
      e4nuConfigurationI( const double EBeam, const unsigned int TargetPdg ) { kEBeam = EBeam ; kTargetPdg = TargetPdg ; }

      virtual ~e4nuConfigurationI() { ; }
      
    protected: 
      bool kUseAllSectors = false ; 
      bool kApplyFiducial = true ; 
      bool kApplyAccWeights = true ; 
      bool kApplyReso = true ; 
      bool kApplyPhiOpeningAngle = false ; 
      bool kUsePhiThetaBand = false ;
      bool kApplyThetaSlice = false; 
      bool kApplyGoodSectorPhiSlice = false ; 
      bool kIsData = false ; 

      // ofset for oscillation studies
      double koffset = 0 ; 

      // Physics 
      double kEBeam = 1.161 ; 
      unsigned int kTargetPdg = 1000060120 ;
    };
  }
}

#endif
