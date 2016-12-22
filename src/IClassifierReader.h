#ifndef ICLASSIFIERREADER_H 
#define ICLASSIFIERREADER_H 1

// Include files

/** @class IClassifierReader IClassifierReader.h
 *  
 *
 *  @author Christian Voss
 *  @date   2014-06-18
 */
class IClassifierReader {
public: 
  /// Standard constructor
  IClassifierReader( )
    : fStatusIsClean( true ) {}
  
  virtual ~IClassifierReader( ) { } ///< Destructor
  
  // return classifier response
  virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;
  
  // returns classifier status
  bool IsStatusClean() const { return fStatusIsClean; }
  
protected:

  bool fStatusIsClean;

private:

};
#endif // ICLASSIFIERREADER_H
