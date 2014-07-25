#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>

#include<string_utils.h>
#include<command_line_options.h>

#include<remesher.h>

std::vector<std::string> tokenize( std::string &in ){
    std::vector<std::string> tokens;
    std::istringstream iss(in);
    std::copy(std::istream_iterator<std::string>(iss),
              std::istream_iterator<std::string>(),
              std::back_inserter<std::vector<std::string> >(tokens));
    return tokens;
}

template< typename real >
bool load_obj( const char *filename, std::vector<real>& coords, std::vector<int>& tris ){
    
    std::ifstream in(filename);
    if( !in.is_open() )
        return false;
    
    std::string line;
    while( std::getline(in,line) ){
        if( line.empty() )
            continue;
        std::istringstream iss(line);
        std::string token;
        
        iss >> token;
        if( token == "v" ){
            real x, y, z;
            iss >> x >> y >> z;
            coords.push_back( x );
            coords.push_back( y );
            coords.push_back( z );
        } else if( token == "f" ){
            std::vector<std::string> tokens = tokenize(line);
            if( tokens.size() != 4 ){
                std::cout << "input not a triangle mesh!" << std::endl;
                return false;
            }
            for( int i=1; i<tokens.size(); i++ ){
                std::replace( tokens[i].begin(), tokens[i].end(), '/', ' ');
                int id = from_str<int>(tokens[i])-1;
                tris.push_back( id );
            }
        }
    }
    
    return true;
}

int main( int argc, const char **argv ){
    char input_file[1024];
    char output_file[1024];
    double edge_len = 1.0/50.0;
    double feat = acos(0.707)*180.0/M_PI;
    int iters=10;
    
    utilities::command_line_options clopts;
    clopts.add_required_parameter( "-input",  ARGUMENT_STRING, input_file,  1, "input file in .obj format" );
    clopts.add_required_parameter( "-output", ARGUMENT_STRING, output_file, 1, "output file in .obj format" );
    clopts.add_optional_parameter( "-size",   ARGUMENT_DOUBLE, &edge_len,   1, "edge length as fraction of bounding box longest side" );
    clopts.add_optional_parameter( "-feat",   ARGUMENT_DOUBLE, &feat,       1, "feature threshold, degrees" );
    clopts.add_optional_parameter( "-iters",  ARGUMENT_INT,    &iters,      1, "number of remeshing iterations to perform");
    if( !clopts.parse( argc, argv ) ){
        return 1;
    }
    
    std::vector<double> coords, rm_coords;
    std::vector<int> tris, rm_tris;
    if( !load_obj( input_file, coords, tris ) )
        return 2;
    
    // compute the longest edge of the
    // input bounding box
    vec3d minim( coords[0], coords[1], coords[2] );
    vec3d maxim = minim;
    for( int i=3; i<coords.size(); i+=3 ){
        vec3d p( coords[i+0], coords[i+1], coords[i+2] );
        minim = minim.min( p );
        maxim = maxim.max( p );
    }
    vec3d delta = maxim-minim;
    double L = std::max( delta[0], std::max( delta[1], delta[2] ) );
    edge_len *= L;
    
    
    remesher_options opts;
    opts["REMESHER_REFINE_FEATURES"]        = "TRUE";
    opts["REMESHER_COARSEN_FEATURES"]       = "TRUE";
    opts["REMESHER_FEATURE_THRESHOLD"]      = to_str(cos(feat*M_PI/180.0));
    opts["REMESHER_MIN_EDGE_LENGTH"]        = to_str(edge_len);
    opts["REMESHER_MAX_EDGE_LENGTH"]        = to_str(edge_len*2.0);
    opts["REMESHER_TARGET_EDGE_LENGTH"]     = to_str(edge_len);
    opts["REMESHER_RELATIVE_EDGE_ERROR"]    = to_str(-1.0);
    opts["REMESHER_ITERATIONS"]             = to_str(iters);
    
    try {
        remesh( opts, coords, tris, rm_coords, rm_tris );
    } catch( const char *err ){
        std::cout << err << std::endl;
        return 1;
    }
    save_obj( output_file, rm_coords, rm_tris );
    
    return 0;
}