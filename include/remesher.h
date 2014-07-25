#ifndef REMESHER_H
#define REMESHER_H

#include<set>
#include<map>
#include<unordered_set>
#include<unordered_map>

#include<vec3.h>
#include<closest_point_on_mesh.h>
#include<string_utils.h>

#include"mesh.h"
using namespace trimesh;

typedef std::pair<vertex*,vertex*> edge;
edge make_edge( vertex* a, vertex *b ){
    return edge( std::min(a,b), std::max(a,b) );
}

void label_features( mesh* m, double thresh ){
    for( mesh::triangle_iterator it=m->triangles_begin(); it!=m->triangles_end(); ++it ){
        triangle* t = *it;
        for( int i=0; i<3; i++ ){
            vertex *a=t->vtx(i);
            vertex *b=t->vtx((i+1)%3);
            triangle* n = m->tri_for_edge(b,a);
            if( !n ){
                a->is_feature() = true;
                b->is_feature() = true;
                m->mark_as_features( a, b );
            } else {
                vec3d n1 = (t->vtx(1)->pos()-t->vtx(0)->pos()).cross(t->vtx(2)->pos()-t->vtx(0)->pos());
                vec3d n2 = (n->vtx(1)->pos()-n->vtx(0)->pos()).cross(n->vtx(2)->pos()-n->vtx(0)->pos());
                n1.normalize();
                n2.normalize();
                if( n1.dot(n2) < thresh ){
                    a->is_feature()=true;
                    b->is_feature()=true;
                    m->mark_as_features( a, b );
                }
            }
        }
    }
}

double min_triangle_angle( vec3d a, vec3d b, vec3d c ){
    vec3d ang = triangle_angles( a, b, c );
    return std::min( ang[0], std::min( ang[1], ang[2] ) );
}

enum {
    FLIP_VALENCE,
    FLIP_MIN_ANGLE,
};

bool edge_should_flip( mesh* m, vertex* a, vertex *b, int mode ){
    triangle* T1 = m->tri_for_edge(a,b);
    triangle* T2 = m->tri_for_edge(b,a);
    if( !T1 || !T2 || m->edge_is_feature(a,b) /*|| a->is_feature() || b->is_feature()*/ )
        return false;
    vertex* c = T1->next(b);
    vertex* d = T2->next(a);
    //if( a->is_feature() || b->is_feature() )
    //    return false;
    
    try {
    vec3d ang1 = triangle_angles( a->pos(), b->pos(), c->pos() );
    vec3d ang2 = triangle_angles( a->pos(), b->pos(), d->pos() );
    double theta_a = ang1[0]+ang2[0];
    double theta_b = ang1[1]+ang2[1];
    if( theta_a >= 0.8*M_PI || theta_b >= 0.8*M_PI )
        return false;
    
    if( mode == FLIP_MIN_ANGLE || a->is_feature() || b->is_feature() ){
        double min_angle1 = std::min( ang1[0], std::min(ang1[1],ang1[2]) );
        double min_angle2 = std::min( ang2[0], std::min(ang2[1],ang2[2]) );
        double min_angle = std::min(min_angle1,min_angle2);
        
        ang1 = triangle_angles( a->pos(), c->pos(), d->pos() );
        ang2 = triangle_angles( b->pos(), d->pos(), c->pos() );
        min_angle1 = std::min( ang1[0], std::min(ang1[1],ang1[2]) );
        min_angle2 = std::min( ang2[0], std::min(ang2[1],ang2[2]) );
        double tmp = std::min(min_angle1,min_angle2);
        return tmp > min_angle;
    }

    // valence flipping strategy
    int va = a->valence()-6;
    int vb = b->valence()-6;
    int vc = c->valence()-6;
    int vd = d->valence()-6;
    int init = va*va + vb*vb + vc*vc + vd*vd;
    va--;
    vb--;
    vc++;
    vd++;
    int final = va*va + vb*vb + vc*vc +vd*vd;
    return final < init;
    } catch( const char *err ){
        return false;
    }
}

void check_closed_manifold( mesh* m ){
    m->check_closed_manifold();
}

void flip_edges( mesh* m, int mode ){
    if( !m->check_mesh() )
        throw "uh oh";
    
    std::vector<edge> edges;
    for( mesh::triangle_iterator it=m->triangles_begin(); it!=m->triangles_end(); ++it ){
        triangle *t = *it;
        edges.push_back( make_edge(t->vtx(0),t->vtx(1)) );
        edges.push_back( make_edge(t->vtx(1),t->vtx(2)) );
        edges.push_back( make_edge(t->vtx(2),t->vtx(0)) );
    }
    for( int i=0; i<edges.size(); i++ ){
        if( edge_should_flip( m, edges[i].first, edges[i].second, mode ) ){
            m->flip_edge( edges[i].first, edges[i].second );
        }
    }
    if( !m->check_mesh() )
        throw "uh oh";

}

template< typename proj >
void split_edges( mesh* m, proj& cp, double minsize, double maxsize, double max_error, bool allow_features ){
    if( !m->check_mesh() )
        throw "uh oh";

    std::vector<edge> edges;
    for( mesh::triangle_iterator it=m->triangles_begin(); it!=m->triangles_end(); ++it ){
        triangle *t = *it;
        edges.push_back( make_edge(t->vtx(0),t->vtx(1)) );
        edges.push_back( make_edge(t->vtx(1),t->vtx(2)) );
        edges.push_back( make_edge(t->vtx(2),t->vtx(0)) );
    }
    for( int i=0; i<edges.size(); i++ ){
        
        vertex* v0 = edges[i].first;
        vertex* v1 = edges[i].second;
        double L = (v1->pos()-v0->pos()).length();
        bool grade_split = false;
        
        bool adapt_split = false;
        
        bool feat_split = !m->edge_is_feature(v0, v1) && (v0->is_feature() && v1->is_feature());
        
        if( max_error > 0.0 ){
            vec3d mid = (v0->pos()+v1->pos())/2.0;
            vec3d n, p = cp( mid, L, n );
            adapt_split = L > minsize && (p-mid).length() > max_error;
        }
        if( L > maxsize || adapt_split || grade_split || feat_split ){
            if( m->edge_is_feature(edges[i].first, edges[i].second) && !allow_features )
                continue;
            m->split_edge(edges[i].first,edges[i].second);
        }
    }
    
    if( !m->check_mesh() )
        throw "uh oh";
}

template< typename proj >
void collapse_edges( mesh* m, proj& cp, double minsize, double maxsize, double max_error, bool allow_features ){
    if( !m->check_mesh() )
        throw "uh oh";

    std::vector<edge> edges;
    for( mesh::triangle_iterator it=m->triangles_begin(); it!=m->triangles_end(); ++it ){
        triangle *t = *it;
        edges.push_back( make_edge(t->vtx(0),t->vtx(1)) );
        edges.push_back( make_edge(t->vtx(1),t->vtx(2)) );
        edges.push_back( make_edge(t->vtx(2),t->vtx(0)) );
    }
    for( int i=0; i<edges.size(); i++ ){
        
        vertex *v0 = edges[i].first;
        vertex *v1 = edges[i].second;
        double L = (v1->pos()-v0->pos()).length();
        
        if( L < minsize || v1->valence() <= 4 ){
            m->collapse_edge(v0,v1,allow_features);
        }
    }

    if( !m->check_mesh() )
        throw "uh oh";
}

double smooth_d = 2.0;

template< typename proj >
void smooth_vertex( mesh* m, proj& prj, vertex* v, bool project ){
    if( v->is_feature() ) return;
    
    vec3d pos(0,0,0);
    vec3d nrm(0.0,0.0,0.0);
    double w = 0.0;
    double min_angle = 180.0;
    for( vertex::triangle_iterator it=v->triangles_begin(); it!=v->triangles_end(); ++it ){
        triangle* t = *it;
        
        min_angle = std::min( min_angle, min_triangle_angle(t->vtx(0)->pos(),t->vtx(1)->pos(),t->vtx(2)->pos()) );
        
        vec3d tnrm = (t->vtx(1)->pos()-t->vtx(0)->pos()).cross(t->vtx(2)->pos()-t->vtx(0)->pos());
        nrm += tnrm;
        double a = tnrm.length();
        for( int j=0; j<3; j++ ){
            
            
            pos += t->vtx(j)->pos();
            w += 1.0;
            
        }
    }
    //pos += 0.001*w*v->pos();
    //w += w*0.001;
    vec3d delta = pos/w - v->pos();
    if( project )
        delta -= nrm*nrm.dot(delta)/(nrm.dot(nrm));

    vec3d oldpos = v->pos();
    
    vec3d n;
    double L = 2.0*delta.length();
    double alpha = 1.0;
    for( int i=0; i<10; i++ ){
        v->pos() += delta*alpha;
        if( project )
            v->pos() = prj( v->pos(), smooth_d, n );
        

        bool valid = true;
        for( vertex::triangle_iterator it=v->triangles_begin(); it!=v->triangles_end(); ++it ){
            triangle *t = *it;
            vertex *v0 = t->next(v);
            vertex *v1 = t->next(v0);
            vec3d oldn = (v0->pos()-oldpos).cross(v1->pos()-oldpos);
            vec3d newn = (v0->pos()-v->pos()).cross(v1->pos()-v->pos());
            if( oldn.dot(newn) <= 0.0 ){
                valid = false;
                break;
            }
        }
        if( valid ){
            break;
        } else {
            v->pos() = oldpos;
            alpha *= 0.75;
        }
        
        
        /*
        double tminangle = 180.0;
        for( vertex::triangle_iterator it=v->triangles_begin(); it!=v->triangles_end(); ++it ){
            triangle* t = *it;
            tminangle = std::min( tminangle, min_triangle_angle( t->vtx(0)->pos(), t->vtx(1)->pos(), t->vtx(2)->pos() ) );
        }
        if( tminangle > min_angle || tminangle > 20.0*M_PI/180 ){
            break;
        } else {
            v->pos() = oldpos;
            alpha *= 0.75;
        }
        */
    }
}

template< typename proj >
void smooth_vertices( mesh* m, proj& prj, bool project ){
    for( mesh::vertex_iterator it=m->vertices_begin(); it!=m->vertices_end(); ++it ){
        vertex *v = *it;
        smooth_vertex(m,prj,v,project);
    }
}

typedef std::map<std::string,std::string> remesher_options;

/*
template< typename T >
T from_str( std::string& tmp ){
    std::istringstream iss(tmp);
    T val;
    iss >> val;
    return val;
}

template<typename T>
std::string to_str( T val ){
    std::ostringstream oss;
    oss << val;
    return oss.str();
}
*/


class remesher {
public:
    remesher( remesher_options opts, const std::vector<double>& icoord, std::vector<int>& itris ) : m_opts(opts), m_icoord(icoord), m_itris(itris) {
        // initialize the remesher
        initialize();
    }
    
    void remesh(){
        for( int i=0; i<m_num_iterations; i++ ){
            std::cout << "iteration [" << i+1 << "/" << m_num_iterations << "]" << std::endl;
            update_edges();
            smooth(true);
            smooth(true);
        }
    }
    
    void get_mesh( std::vector<double>& coords, std::vector<int>& tris ){
        m_mesh->get_mesh( coords, tris );
    }
    
private:
    
    double target_size( vec3d pos ){
        return m_target_edge_length;
    }
    
    bool update_mesh_guarded( std::vector<triangle*>& rmtri, std::vector<vertex*>& addtrivtx ){
        // remove the triangles from the mesh, but store the
        // vertex indices that define them
        std::vector<vertex*> rmtrivtx;
        for( int i=0; i<rmtri.size(); i++ ){
            rmtrivtx.push_back( rmtri[i]->vtx(0) );
            rmtrivtx.push_back( rmtri[i]->vtx(1) );
            rmtrivtx.push_back( rmtri[i]->vtx(2) );
            m_mesh->remove_triangle( rmtri[i] );
        }
   
        // now attempt to add the triangles, one by one,
        // checking to see if they are successfully added
        triangle *t;
        bool failed = false;
        std::vector<triangle*> addtri;
        for( int i=0; i<addtrivtx.size(); i+=3 ){
            t = m_mesh->add_triangle( addtrivtx[i+0], addtrivtx[i+1], addtrivtx[i+2] );
            if( !t ){
                failed=true;
                break;
            } else {
                addtri.push_back(t);
            }
        }
        if( failed ){
            for( int i=0; i<addtri.size(); i++ ){
                m_mesh->remove_triangle( addtri[i] );
            }
            for( int i=0; i<rmtrivtx.size(); i+=3 ){
                m_mesh->add_triangle( rmtrivtx[i+0], rmtrivtx[i+1], rmtrivtx[i+2] );
            }
            return false;
        }
        return true;
    }
    
    void smooth(bool proj){
        for( mesh::vertex_iterator it=m_mesh->vertices_begin(); it!=m_mesh->vertices_end(); ++it ){
            vertex* v = *it;
            if( v->is_feature() )
                continue;
            vec3d vpos(0.0,0.0,0.0);
            double wgt=0.0;
            for( vertex::triangle_iterator tit=v->triangles_begin(); tit!=v->triangles_end(); ++tit ){
                triangle* t = *tit;
                double w = (t->vtx(1)->pos()-t->vtx(0)->pos()).cross(t->vtx(2)->pos()-t->vtx(0)->pos()).length();
                vpos += w*(t->vtx(0)->pos()+t->vtx(1)->pos()+t->vtx(2)->pos())/3.0;
                wgt += w;
            }
            vec3d nrm;
            vpos /= wgt;
            if( proj )
                vpos = m_proj(vpos,m_closest_point_search_radius,nrm);
            
            // project vpos to the surface
            bool do_it=true;
            for( vertex::triangle_iterator tit=v->triangles_begin(); tit!=v->triangles_end(); ++tit ){
                triangle* t = *tit;
                vertex *a = v;
                vertex *b = t->next(v);
                vertex *c = t->next(b);
                vec3d curnrm = (b->pos()-a->pos()).cross(c->pos()-a->pos());
                vec3d newnrm = (b->pos()-vpos).cross(c->pos()-vpos);
                if( curnrm.dot(newnrm) > 0.0 ){
                    double oldang = min_triangle_angle( a->pos(), b->pos(), c->pos() );
                    double newang = min_triangle_angle( vpos, b->pos(), c->pos() );
                    if( newang < oldang ){
                        //do_it = false;
                        //break;
                    }
                }
            }
            if( do_it ){
                v->pos() = vpos;
            }
        }
    }
    
    void update_edges(){
        int e[][2] = { {0,1},{1,2},{2,0} };
        // build a list of edges
        /*
        std::set< std::pair<double,edge> > edge_set;
        for( mesh::triangle_iterator it=m_mesh->triangles_begin(); it!=m_mesh->triangles_end(); ++it ){
            triangle* tri = *it;
            for( int j=0; j<3; j++ ){
                if( tri->vtx(e[j][0]) < tri->vtx(e[j][1]) ){
                    edge_set.insert( std::pair<double,edge>((tri->vtx(e[j][0])->pos()-tri->vtx(e[j][1])->pos()).length(),edge(tri->vtx(e[j][0]),tri->vtx(e[j][1]) ) ) );
                }
            }
        }
        */
        std::vector< edge > edge_list;
        //for( std::set< std::pair<double,edge> >::iterator it=edge_set.begin(); it!=edge_set.end(); ++it ){
        //    edge_list.push_back( (*it).second );
        //}
        for( mesh::triangle_iterator it=m_mesh->triangles_begin(); it!=m_mesh->triangles_end(); ++it ){
            triangle* tri = *it;
            for( int j=0; j<3; j++ ){
                if( tri->vtx(e[j][0]) < tri->vtx(e[j][1]) ){
                    edge_list.push_back( edge(tri->vtx(e[j][0]),tri->vtx(e[j][1]) ) );
                }
            }
        }
        
        // now loop over the edges and decide what to do
        for( int i=0; i<edge_list.size(); i++ ){
            
            // first grab the edge vertices and check that
            // they are still in the mesh
            vertex* a = edge_list[i].first;
            vertex* b = edge_list[i].second;
            if( !m_mesh->has_vertex(a) || !m_mesh->has_vertex(b) )
                continue;
            
            // now get the adjacent triangles and check
            // that they are also still in the mesh
            triangle* t0 = m_mesh->tri_for_edge(a,b);
            triangle* t1 = m_mesh->tri_for_edge(b,a);
            if( !m_mesh->has_triangle(t0) || !m_mesh->has_triangle(t1) )
                continue;
            
            // now get the opposing vertices for the
            // triangles in order to check if a flip
            // should be performed
            vertex *c = t0->other(a,b);
            vertex *d = t1->other(a,b);
            
            // compute the angles in the triangle to
            // make sure that the unfolded quad formed
            // by t0 and t1 is convex
            vec3d t0_angles = triangle_angles( a->pos(), b->pos(), c->pos() );
            vec3d t1_angles = triangle_angles( a->pos(), b->pos(), d->pos() );
            
            // check if the edge is a feature edge
            bool is_feature_edge = m_mesh->edge_is_feature(a,b);
            
            // check for convexity, and if so, consider
            // flipping the edge ab -> cd
            double convex_threshold = 0.9;
            if( t0_angles[0]+t1_angles[0] < convex_threshold*M_PI && t0_angles[1]+t1_angles[1] < convex_threshold*M_PI ){
                // sum of angles at vertex a and b are less than
                // 180 degrees, so [a,b] is a flip candidate
                
                int va = a->valence();
                int vb = b->valence();
                int vc = c->valence();
                int vd = d->valence();
                int init_L2 = va*va + vb*vb + vc*vc + vd*vd;
                va--;
                vb--;
                vc++;
                vd++;
                int final_L2 = va*va + vb*vb + vc*vc + vd*vd;
                
                vec3d newangles1 = triangle_angles( a->pos(), d->pos(), c->pos() );
                vec3d newangles2 = triangle_angles( b->pos(), c->pos(), d->pos() );
                newangles1 = newangles1.min(newangles2);
                double newminangle = std::min( newangles1[0], std::min(newangles1[1],newangles1[2]));
                newangles1 = t0_angles.min(t1_angles);
                double oldminangle = std::min( newangles1[0], std::min(newangles1[2],newangles1[2]));
                
                // check if the edge should flip, default condition is that
                // the edge is not a feature edge
                bool edge_should_flip = !is_feature_edge;
                
                // minimum angle should improve as a result of an edge flip
                edge_should_flip &= init_L2 > final_L2 && newminangle > 5*M_PI/180.0;
                if( a->is_feature() || b->is_feature() )
                    edge_should_flip &= newminangle > oldminangle /*&& newminangle > 5*M_PI/180.0*/;
                
                if( edge_should_flip ){
                    std::vector<vertex*>   addtrivtx;
                    std::vector<triangle*> rmtri;
                    // remove the two triangles
                    rmtri.push_back( t0 );
                    rmtri.push_back( t1 );
                    
                    // add the flipped triangles
                    addtrivtx.push_back(a); addtrivtx.push_back(d); addtrivtx.push_back(c);
                    addtrivtx.push_back(b); addtrivtx.push_back(c); addtrivtx.push_back(d);
                    
                    // try to update the mesh
                    if( !update_mesh_guarded( rmtri, addtrivtx ) ){
                    } else {
                    }
                    
                    // if the edge flipped then don't bother
                    // trying the rest of the options for
                    // the edge
                    continue;
                }
            }
            
            // compute the squared length of the edge
            // and its midpoint in preparation for
            // deciding whether to collapse or
            // subdivide the edge, also retrieve the
            // target size from the sizing field.
            double L  = (a->pos()-b->pos()).length_squared();
            vec3d mid = (a->pos()+b->pos())/2.0;
            double ts = target_size(mid);
            
            double max_edge_length_sq = (2.0*ts)*(2.0*ts);
            if( L > max_edge_length_sq ){
                // edge is longer than the maximum edge length, see if it
                // can be subdivided
                bool edge_should_subdivide = true;
                
                // do not subdivide feature edges unless explicitly permitted
                edge_should_subdivide &= !is_feature_edge || (is_feature_edge && m_split_features);
                
                // now try to subdivide the
                if( edge_should_subdivide ){
                    
                    // generate the new midpoint vertex
                    vec3d nrm;
                    vertex *m = m_mesh->add_vertex( (a->pos()+b->pos())*0.5 );
                    m->pos() = m_proj( m->pos(), m_closest_point_search_radius, nrm );
                    
                    // plan the mesh modifications
                    std::vector<triangle*> rmtri;
                    std::vector<vertex*> addtrivtx;
                    rmtri.push_back( t0 );
                    rmtri.push_back( t1 );
                    addtrivtx.push_back( a ); addtrivtx.push_back( m ); addtrivtx.push_back( c );
                    addtrivtx.push_back( a ); addtrivtx.push_back( d ); addtrivtx.push_back( m );
                    addtrivtx.push_back( b ); addtrivtx.push_back( m ); addtrivtx.push_back( d );
                    addtrivtx.push_back( b ); addtrivtx.push_back( c ); addtrivtx.push_back( m );
                    
                    vec3d angles = triangle_angles( a->pos(), m->pos(), c->pos() );
                    angles = angles.min( triangle_angles( a->pos(), d->pos(), m->pos() ) );
                    angles = angles.min( triangle_angles( b->pos(), m->pos(), d->pos() ) );
                    angles = angles.min( triangle_angles( b->pos(), c->pos(), m->pos() ) );
                    double min_angle = std::min(angles[0],std::min(angles[1],angles[2]));
                    
                    // perform the mesh update and delete
                    // the new vertex if it fails
                    if( /*min_angle > 1.0*M_PI/180.0 &&*/ update_mesh_guarded( rmtri, addtrivtx ) ){
                        if( is_feature_edge ){
                            m_mesh->mark_as_features( a, m );
                            m_mesh->mark_as_features( m, b );
                        }
                    } else {
                        m_mesh->remove_vertex(m);
                    }
                    
                    continue;
                }
            }
            
            // check if the edge is too short, in which case
            // it should be collapsed, pending validity checks
            double min_edge_length_sq = (0.5*ts)*(0.5*ts);
            if( L < min_edge_length_sq ){
                // edge is shorter than the minimum edge length, see
                // if it should be collapsed.
                bool edge_should_collapse = true;
                
                // check if edge is a feature and collapsible, or not a feature
                edge_should_collapse &= !is_feature_edge || (is_feature_edge && m_coarsen_features);
                
                if( !is_feature_edge && b->is_feature() )
                    edge_should_collapse = false;
                
                if( a->valence() <= 3 || b->valence() <= 3 )
                    edge_should_collapse &= true;
                
                if( a->valence() == 5 && b->valence() == 5 )
                    edge_should_collapse &= true;
                
                if( edge_should_collapse ){
                    bool would_flip_tri = false;
                    std::vector<triangle*> rmtri;
                    std::vector<vertex*> addtrivtx;
                    std::vector<vertex*> newfeature;
                    for( vertex::triangle_iterator it=b->triangles_begin(); it!=b->triangles_end(); ++it ){
                        triangle* t = *it;
                        // if the triangle uses vertex a, it will be
                        // deleted in the edge collapse, add it to
                        // the removed triangles and continue
                        if( t->has_vertex(a) ){
                            rmtri.push_back( t );
                            continue;
                        }
                        // triangle touches vertex b only, check the
                        // normal before and after the collapse and
                        // make sure the triangle would not flip, also
                        // check that the minimum angle is larger than
                        // 10 degrees, since we don't want to generate
                        // degenerate triangles
                        vertex* t0 = b;
                        vertex* t1 = t->next(b);
                        vertex* t2 = t->next(t1);
                        vec3d curnrm = (t1->pos()-t0->pos()).cross(t2->pos()-t0->pos());
                        vec3d newnrm = (t1->pos()-a->pos()).cross(t2->pos()-a->pos());
                        vec3d newang = triangle_angles( a->pos(), t1->pos(), t2->pos() );
                        double minang = std::min(newang[0],std::min(newang[1],newang[2]));
                        if( newnrm.dot(curnrm) <= 0.0 /*|| minang < 1.0*M_PI/180.0*/ ){
                            would_flip_tri = true;
                            break;
                        }
                        
                        if( m_mesh->edge_is_feature(b,t1) ){
                            newfeature.push_back(a);
                            newfeature.push_back(t1);
                        }
                        if( m_mesh->edge_is_feature(b,t2) ){
                            newfeature.push_back(a);
                            newfeature.push_back(t2);
                        }
                        
                        // add the triangle to the removed list and
                        // add the new triangle vertices to the
                        // new triangles vertex list
                        rmtri.push_back( t );
                        addtrivtx.push_back( a );
                        addtrivtx.push_back( t1 );
                        addtrivtx.push_back( t2 );
                    }
                    
                    if(!would_flip_tri){
                        if( update_mesh_guarded( rmtri, addtrivtx) ){
                            m_mesh->remove_vertex(b);
                            for( int j=0; j<newfeature.size(); j+=2 ){
                                m_mesh->mark_as_features( newfeature[j+0], newfeature[j+1] );
                            }
                        }
                    }
                    
                    continue;
                }
            }
        }
    }
    
    void label_features(){
        double thresh = cos(m_feature_threshold);
        for( mesh::triangle_iterator it=m_mesh->triangles_begin(); it!=m_mesh->triangles_end(); ++it ){
            triangle* t = *it;
            for( int i=0; i<3; i++ ){
                vertex *a=t->vtx(i);
                vertex *b=t->vtx((i+1)%3);
                triangle* n = m_mesh->tri_for_edge(b,a);
                if( !n ){
                    a->is_feature() = true;
                    b->is_feature() = true;
                    m_mesh->mark_as_features( a, b );
                } else {
                    vec3d n1 = (t->vtx(1)->pos()-t->vtx(0)->pos()).cross(t->vtx(2)->pos()-t->vtx(0)->pos());
                    vec3d n2 = (n->vtx(1)->pos()-n->vtx(0)->pos()).cross(n->vtx(2)->pos()-n->vtx(0)->pos());
                    n1.normalize();
                    n2.normalize();
                    if( n1.dot(n2) < thresh ){
                        a->is_feature()=true;
                        b->is_feature()=true;
                        m_mesh->mark_as_features( a, b );
                    }
                }
            }
        }
    }
    
    void initialize(){
        // handle the input options and perform checking
        if( m_opts["REMESHER_FEATURE_THRESHOLD"] == "" )   m_opts["REMESHER_FEATURE_THRESHOLD"]   = to_str<double>(0.5);
        if( m_opts["REMESHER_COARSEN_FEATURES"]  == "" )   m_opts["REMESHER_COARSEN_FEATURES"]    = "TRUE";
        if( m_opts["REMESHER_REFINE_FEATURES"]   == "" )   m_opts["REMESHER_REFINE_FEATURES"]     = "TRUE";
        if( m_opts["REMESHER_TARGET_EDGE_LENGTH"]  == "" ) m_opts["REMESHER_TARGET_EDGE_LENGTH"]  = to_str<double>(1.0);
        if( m_opts["REMESHER_ITERATIONS"] == "" )          m_opts["REMESHER_ITERATIONS"]          = to_str<int>(10);
        
        m_feature_threshold  = from_str<double>(m_opts["REMESHER_FEATURE_THRESHOLD"]);
        m_target_edge_length = from_str<double>(m_opts["REMESHER_TARGET_EDGE_LENGTH"]);
        m_num_iterations     = from_str<int>(m_opts["REMESHER_ITERATIONS"]);
        m_split_features     = m_opts["REMESHER_REFINE_FEATURES"] == "TRUE";
        m_coarsen_features   = m_opts["REMESHER_COARSEN_FEATURES"] == "TRUE";
        
        // create the mesh and store a list of vertex pointers
        // and positions to build the triangles and closest point
        // search structure
        m_mesh = new mesh();
        std::vector<vertex*> vtx;
        std::vector<vec3d>   vtxpos;
        for( int i=0; i<m_icoord.size(); i+=3 ){
            vec3d pos( m_icoord[i+0], m_icoord[i+1], m_icoord[i+2] );
            vtxpos.push_back( pos );
            vtx.push_back( m_mesh->add_vertex( pos ) );
            if( i == 0 ){
                m_mesh_minim = pos;
                m_mesh_maxim = pos;
            } else {
                m_mesh_minim = m_mesh_minim.min( pos );
                m_mesh_maxim = m_mesh_maxim.max( pos );
            }
        }
        m_proj.initialize( vtxpos, m_itris );
        m_closest_point_search_radius = 2.0*m_target_edge_length;
        
        // now add the mesh triangles
        for( int i=0; i<m_itris.size(); i+=3 ){
            m_mesh->add_triangle( vtx[m_itris[i+0]], vtx[m_itris[i+1]], vtx[m_itris[i+2]] );
        }
        label_features();
    }
    
    remesher_options                            m_opts;
    const std::vector<double>&                  m_icoord;
    const std::vector<int>&                     m_itris;
    geom::closest_point_on_mesh<double,vec3d>   m_proj;
    
    double m_feature_threshold;
    double m_target_edge_length;
    int    m_num_iterations;
    bool   m_split_features;
    bool   m_coarsen_features;
    
    mesh*  m_mesh;
    vec3d  m_mesh_minim;
    vec3d  m_mesh_maxim;
    double m_closest_point_search_radius;
};


void remesh( remesher_options opts, const std::vector<double>& icoord, std::vector<int>& itris, std::vector<double>&ocoord, std::vector<int>& otris ){
    // set up the input options, setting defaults if none are specified
    if( opts["REMESHER_FEATURE_THRESHOLD"] == "" )   opts["REMESHER_FEATURE_THRESHOLD"]   = to_str<double>(0.5);
    if( opts["REMESHER_COARSEN_FEATURES"]  == "" )   opts["REMESHER_COARSEN_FEATURES"]    = "TRUE";
    if( opts["REMESHER_REFINE_FEATURES"]   == "" )   opts["REMESHER_REFINE_FEATURES"]     = "TRUE";
    if( opts["REMESHER_MIN_EDGE_LENGTH"]   == "" )   opts["REMESHER_MIN_EDGE_LENGTH"]     = to_str<double>(1.0);
    if( opts["REMESHER_MAX_EDGE_LENGTH"]   == "" )   opts["REMESHER_MAX_EDGE_LENGTH"]     = to_str<double>(2.0);
    if( opts["REMESHER_RELATIVE_EDGE_ERROR"] == "" ) opts["REMESHER_RELATIVE_EDGE_ERROR"] = to_str<double>(-1000.0);
    if( opts["REMESHER_ITERATIONS"] == "" )          opts["REMESHER_ITERATIONS"]          = to_str<int>(10);
    double feature_thresh   = from_str<double>(opts["REMESHER_FEATURE_THRESHOLD"]);
    double min_edge_length  = from_str<double>(opts["REMESHER_MIN_EDGE_LENGTH"]);
    double max_edge_length  = from_str<double>(opts["REMESHER_MAX_EDGE_LENGTH"]);
    double rel_edge_error   = from_str<double>(opts["REMESHER_RELATIVE_EDGE_ERROR"]);
    int    iterations       = from_str<int>(opts["REMESHER_ITERATIONS"]);
    bool   split_features   = opts["REMESHER_REFINE_FEATURES"] == "TRUE";
    bool   coarsen_features = opts["REMESHER_COARSEN_FEATURES"] == "TRUE";

    mesh m;
    vec3d minim( icoord[0], icoord[1], icoord[2] );
    vec3d maxim( icoord[0], icoord[1], icoord[2] );
    std::vector<vertex*> v;
    std::vector<vec3d> tv;
    for( int i=0; i<icoord.size(); i+=3 ){
        v.push_back( m.add_vertex( vec3d(icoord[i+0], icoord[i+1],icoord[i+2]) ) );
        tv.push_back( vec3d(icoord[i+0],icoord[i+1],icoord[i+2]));
        minim = minim.min( vec3d(icoord[i+0], icoord[i+1], icoord[i+2]) );
        maxim = maxim.max( vec3d(icoord[i+0], icoord[i+1], icoord[i+2]) );

    }
    smooth_d = max_edge_length*2.0;
    for( int i=0; i<itris.size(); i+=3 ){
        m.add_triangle( v[itris[i+0]], v[itris[i+1]], v[itris[i+2]] );
    }
    
    remesher rm( opts, icoord, itris );
    rm.remesh();
    rm.get_mesh( ocoord, otris );
    return;
    
    
    geom::closest_point_on_mesh<double,vec3d> cp;
    cp.initialize( tv, itris );
    
    label_features( &m, feature_thresh);

    for( int i=0; i<10; i++ ){
        flip_edges( &m, FLIP_MIN_ANGLE );
    }
    
    for( int i=0; i<iterations; i++ ){
        std::cout << "remesh iteration [" << i+1 << "/" << iterations << "]" << std::endl;
        flip_edges( &m, FLIP_VALENCE );
        split_edges( &m, cp, min_edge_length, max_edge_length, -1.0, split_features );
        smooth_vertices( &m, cp, true );
        flip_edges( &m, FLIP_VALENCE );
        smooth_vertices( &m, cp, true );
        flip_edges( &m, FLIP_MIN_ANGLE );
        collapse_edges(&m,cp,min_edge_length,max_edge_length,rel_edge_error,coarsen_features);
    }
    flip_edges( &m, FLIP_MIN_ANGLE );
    smooth_vertices( &m, cp, true );
    
    /*
    int subdiv_iters=5;
    for( int i=0; i<subdiv_iters; i++ ){
        std::cout << "split iteration [" << i << "/" << subdiv_iters << "]" << std::endl;
        split_edges( &m, cp, min_edge_length, max_edge_length, rel_edge_error, split_features );
        flip_edges( &m, FLIP_VALENCE );
        //collapse_edges(&m,cp,min_edge_length,max_edge_length,rel_edge_error,coarsen_features);
        //flip_edges( &m, FLIP_VALENCE );
        smooth_vertices( &m, cp, true );
    }
    
    int remesh_iters = 10;
    for( int i=0; i<remesh_iters; i++ ){
        std::cout << "remesh iteration [" << i << "/" << remesh_iters << "]" << std::endl;
        smooth_vertices( &m, cp, true );
        collapse_edges(&m,cp,min_edge_length,max_edge_length,rel_edge_error,coarsen_features);
        flip_edges(&m, FLIP_VALENCE);
        smooth_vertices( &m, cp, true );
        //split_edges( &m, cp, min_edge_length, max_edge_length, rel_edge_error, split_features );
        flip_edges( &m, FLIP_MIN_ANGLE );
    }
    */

    ocoord.clear();
    otris.clear();
    m.get_mesh( ocoord, otris );

}

#endif
