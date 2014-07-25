#ifndef MESH_H
#define MESH_H

#include<set>
#include<map>
#include<vector>
#include<unordered_set>
#include<fstream>
#include<iostream>

#include<vec3.h>

namespace trimesh {
  
    template< typename real >
    bool save_obj( const char *filename, std::vector<real>& coords, std::vector<int>& tris ){
        std::ofstream out(filename);
        for( int i=0; i<coords.size(); i+=3 ){
            out << "v " << coords[i+0] << " " << coords[i+1] << " " << coords[i+2] << std::endl;
        }
        for( int i=0; i<tris.size(); i+=3 ){
            out << "f " << tris[i+0]+1 << " " << tris[i+1]+1 << " " << tris[i+2]+1 << std::endl;
        }
        out.close();
        return true;
    }
    
    typedef geom::vec3<double> vec3d;
    
    class mesh;
    class vertex;
    class triangle;
    
    vec3d triangle_angles( vec3d A, vec3d B, vec3d C ){
        double a = (B-C).length();
        double b = (A-C).length();
        double c = (A-B).length();
        double cosalpha = std::max( -1.0, std::min( 1.0, (b*b+c*c-a*a)/(2.0*b*c) ) );
        double cosbeta  = std::max( -1.0, std::min( 1.0, (a*a+c*c-b*b)/(2.0*a*c) ) );
        double cosgamma = std::max( -1.0, std::min( 1.0, (a*a+b*b-c*c)/(2.0*a*b) ) );
        vec3d ret( acos(cosalpha),acos(cosbeta),acos(cosgamma));
        if( fabs(ret[0]+ret[1]+ret[2]-M_PI) > 1e-5 ){
            std::cout << a << ", " << b << ", " << c << ", " << cosalpha << ", " << cosbeta << ", " << cosgamma << ", " << ret[0]+ret[1]+ret[2]-M_PI << std::endl;
            //throw "triangle angle formula error";
            return vec3d(0,0,0);
        }
        return ret;
    }
    
    class vertex {
    public:
        typedef std::set<triangle*>::iterator triangle_iterator;
        
        vertex( mesh *m, vec3d p ) : m_mesh(m), m_pos(p) {
            m_is_feature=false;
        }
        
        bool& is_feature(){ return m_is_feature; }
        
        int &feature_edge_count(){ return m_feat_count; }
        
        void add_triangle( triangle* t ){
            m_tris.insert( t );
        }
        
        void remove_triangle( triangle* t ){
            m_tris.erase(t);
        }
        
        bool has_triangle( triangle* t ){
            return m_tris.find(t) != m_tris.end();
        }
        
        int valence(){
            return m_tris.size();
        }
        
        triangle_iterator triangles_begin(){
            return m_tris.begin();
        }
        triangle_iterator triangles_end(){
            return m_tris.end();
        }
        
        vec3d& pos(){
            return m_pos;
        }
    private:
        mesh*               m_mesh;
        vec3d               m_pos;
        std::set<triangle*> m_tris;
        bool m_is_feature;
        int m_feat_count;
    };
    
    class triangle {
    public:
        triangle( mesh *m, vertex* a, vertex *b, vertex* c ) : m_mesh(m) {
            m_vtx[0]=a; m_vtx[1]=b; m_vtx[2]=c;
        }
        
        virtual ~triangle(){
        }
        
        bool has_vertex( vertex* v ){
            return m_vtx[0]==v || m_vtx[1] == v || m_vtx[2] == v;
        }
        
        vertex* vtx( int id ){
            return m_vtx[id];
        }
        
        vertex* next( vertex* v ){
            if( m_vtx[0] == v ) return m_vtx[1];
            if( m_vtx[1] == v ) return m_vtx[2];
            if( m_vtx[2] == v ) return m_vtx[0];
            return NULL;
        }
        
        vertex* other( vertex* a, vertex* b ){
            int cnt=0;
            vertex* oth=NULL;
            if( m_vtx[0] == a || m_vtx[0] == b ) cnt++;
            else oth=m_vtx[0];
            if( m_vtx[1] == a || m_vtx[1] == b ) cnt++;
            else oth=m_vtx[1];
            if( m_vtx[2] == a || m_vtx[2] == b ) cnt++;
            else oth=m_vtx[2];
            return cnt == 2 ? oth : NULL;
        }
        
        bool shares_edge( triangle* t ){
            return (t->has_vertex( vtx(0) ) && t->has_vertex( vtx(1) ))
                || (t->has_vertex( vtx(1) ) && t->has_vertex( vtx(2) ))
                || (t->has_vertex( vtx(2) ) && t->has_vertex( vtx(0) ));
        }
    private:
        mesh*   m_mesh;
        vertex* m_vtx[3];
    };
    
    class mesh {
    public:
        typedef std::unordered_set<vertex*> vertex_set;
        typedef std::unordered_set<triangle*> triangle_set;
        typedef typename vertex_set::iterator vertex_iterator;
        typedef typename triangle_set::iterator triangle_iterator;
        
        vertex* add_vertex( vec3d pos ){
            vertex* v = new vertex(this,pos);
            m_verts.insert(v);
            return v;
        }
        
        void remove_vertex( vertex *v ){
            if( v->valence() != 0 ){
                fail_and_die("tried to remove a connected vertex");
            }
            m_verts.erase(v);
            delete v;
        }
        
        bool has_vertex( vertex* v ){
            return m_verts.find(v) != m_verts.end();
        }
        
        triangle* add_triangle( vertex* a, vertex *b, vertex *c ){
            if( tri_for_edge(a,b) || tri_for_edge(b,c) || tri_for_edge(c,a) ){
                //fail_and_die("add_triangle: new triangle would create non manifold edge");
                return NULL;
            }
            if( a == b || a == c || b == c )
                fail_and_die("add_triangle: triangle has duplicated vertices");
            
            if( !has_vertex(a) || !has_vertex(b) || !has_vertex(c) )
                fail_and_die("add_triangle: triangle references unknown vertices" );
            
            for( vertex::triangle_iterator it=a->triangles_begin(); it!=a->triangles_end(); ++it ){
                triangle *t=*it;
                if( t->has_vertex(b) && t->has_vertex(c) )
                    fail_and_die("add_triangle: triangle would be duplicated" );
            }
            for( vertex::triangle_iterator it=b->triangles_begin(); it!=b->triangles_end(); ++it ){
                triangle *t=*it;
                if( t->has_vertex(a) && t->has_vertex(c) )
                    fail_and_die("add_triangle: triangle would be duplicated" );
            }
            for( vertex::triangle_iterator it=c->triangles_begin(); it!=c->triangles_end(); ++it ){
                triangle *t=*it;
                if( t->has_vertex(a) && t->has_vertex(b) )
                    fail_and_die("add_triangle: triangle would be duplicated" );
            }
            
            triangle* t = new triangle(this,a,b,c);
            m_tris.insert(t);
            a->add_triangle(t);
            b->add_triangle(t);
            c->add_triangle(t);
            m_tri_for_edge[ halfedge(a,b) ] = t;
            m_tri_for_edge[ halfedge(b,c) ] = t;
            m_tri_for_edge[ halfedge(c,a) ] = t;
            if( !check_mesh_pointers(t) )
                fail_and_die("mesh pointers incorrect");
            
            return t;
        }
        
        void remove_triangle( triangle* t ){
            vertex *a,*b,*c;
            a = t->vtx(0);
            b = t->vtx(1);
            c = t->vtx(2);
            
            m_tris.erase( t );
            a->remove_triangle(t);
            b->remove_triangle(t);
            c->remove_triangle(t);
            m_tri_for_edge.erase( halfedge(a,b) );
            m_tri_for_edge.erase( halfedge(b,c) );
            m_tri_for_edge.erase( halfedge(c,a) );
            delete t;
            if( tri_for_edge(a,b) || tri_for_edge(b,c) || tri_for_edge(c,a) )
                fail_and_die("remove_triangle: local error in mesh");
        }
        
        bool has_triangle( triangle* t ){
            return m_tris.find(t) != m_tris.end();
        }
        
        vertex_iterator vertices_begin(){
            return m_verts.begin();
        }
        vertex_iterator vertices_end(){
            return m_verts.end();
        }
        triangle_iterator triangles_begin(){
            return m_tris.begin();
        }
        triangle_iterator triangles_end(){
            return m_tris.end();
        }
        
        triangle* tri_for_edge( vertex* a, vertex* b ){
            std::map< halfedge, triangle*>::iterator it=m_tri_for_edge.find( halfedge(a,b) );
            if( it != m_tri_for_edge.end() ){
                triangle* t = it->second;
                if( m_tris.find(t) != m_tris.end() )
                    return t;
            }
            return NULL;
        }
        
        void flip_edge( vertex* a, vertex* b ){
            if( a->valence() == 3 || b->valence() == 3 )
                return;
            
            if( !check_local_mesh(a) || !check_local_mesh(b)  )
                fail_and_die("flip_edge: local error in mesh");
            
            triangle* T1 = tri_for_edge(a,b);
            triangle* T2 = tri_for_edge(b,a);
            if( !T1 || !T2 ) return;
            vertex* c = T1->other(a,b);
            vertex* d = T2->other(a,b);
            if( !check_local_mesh(c) || !check_local_mesh(d) )
                fail_and_die("flip_edge: local error in mesh");
            
            remove_triangle( T1 );
            remove_triangle( T2 );
            T1 = add_triangle( c, a, d );
            T2 = add_triangle( d, b, c );
            if( !T1 || !T2 ){
                if( T1 ) remove_triangle(T1);
                if( T2 ) remove_triangle(T2);
                add_triangle( a, b, c );
                add_triangle( b, a, d );
            }
            
            if( !check_local_mesh(a) || !check_local_mesh(b) || !check_local_mesh(c) || !check_local_mesh(d) )
                fail_and_die("flip_edge: local error in mesh");
        }
        
        void extract_triangles_near_edge( vertex* a, vertex* b, std::vector<double>& coord, std::vector<int>& tris ){
            int cnt=0;
            for( vertex::triangle_iterator it=a->triangles_begin(); it!=a->triangles_end(); ++it ){
                triangle *t = *it;
                for( int i=0; i<3; i++ ){
                    vertex* v = t->vtx(i);
                    coord.push_back( v->pos()[0] );
                    coord.push_back( v->pos()[1] );
                    coord.push_back( v->pos()[2] );
                    tris.push_back(cnt++);
                }
            }
            for( vertex::triangle_iterator it=b->triangles_begin(); it!=b->triangles_end(); ++it ){
                triangle* t = *it;
                if( !t->has_vertex(a) ){
                    for( int i=0; i<3; i++ ){
                        vertex* v = t->vtx(i);
                        coord.push_back( v->pos()[0] );
                        coord.push_back( v->pos()[1] );
                        coord.push_back( v->pos()[2] );
                        tris.push_back(cnt++);
                    }
                }
            }
        }
        
        void collapse_edge( vertex* a, vertex* b, bool allow_features ){
            // do not collapse if the edge is a feature edge
            if( edge_is_feature(a, b) && !allow_features )
                return;
            // if the edge is not a feature edge, but b is a feature
            // vertex, do not collapse because this will change
            // features associated with b
            if( !edge_is_feature(a,b) && (b->is_feature() || a->is_feature()) )
                return;

            bool feat_edge = edge_is_feature(a,b);
            int  feat_edge_count=0;
            
            
            triangle* T1 = tri_for_edge(a,b);
            triangle* T2 = tri_for_edge(b,a);
            if( !T1 || !T2 )
                return;
            
            // needed?
            //if( a->valence()-2+b->valence()-2 <= 3 || b->valence() == 3 || T1->other(a,b)->valence() == 3 || T2->other(a,b)->valence() == 3 )
            //    return;
   
            vec3d newpos = a->pos();
            if( (feat_edge && feat_edge_count < 6) || (!a->is_feature() && !b->is_feature()) )
                newpos = (a->pos()+b->pos())/2.0;
            
            
            std::vector<triangle*> remtri;
            std::vector<triangle*> newtri;
            std::vector<vertex*> newtrivtx;
            std::vector<vertex*> oldtrivtx;
            vertex* v = a; //add_vertex( ( a->pos() + b->pos() )/2.0 );
            for( vertex::triangle_iterator it=a->triangles_begin(); it!=a->triangles_end(); ++it ){
                triangle* t = *it;
                vertex* v1 = t->next(a);
                vertex* v2 = t->next(v1);
                
                feat_edge_count += edge_is_feature(a,v1);
                feat_edge_count += edge_is_feature(a,v2);
                
                vec3d nstart = (v1->pos()-a->pos()).cross(v2->pos()-a->pos());
                vec3d nend   = (v1->pos()-newpos).cross(v2->pos()-newpos);
                if( nstart.dot(nend) <= 0.0 )
                    return;
                
                if( !t->has_vertex(b) ){
                    if( v1 == a || v1 == b )
                        fail_and_die("uhoh");
                    if( v2 == a || v2 == b || v2 == v1 )
                        fail_and_die("uhoh");
                    
                    //newtrivtx.push_back(v);
                    //newtrivtx.push_back(v1);
                    //newtrivtx.push_back(v2);
                } else {
                    oldtrivtx.push_back(a);
                    oldtrivtx.push_back(v1);
                    oldtrivtx.push_back(v2);
                    remtri.push_back(t);
                }
            }
            std::vector< halfedge > newfeat;
            for( vertex::triangle_iterator it=b->triangles_begin(); it!=b->triangles_end(); ++it ){
                triangle *t = *it;
                if( !t->has_vertex(a) ){
                    vertex* v1 = t->next(b);
                    vertex* v2 = t->next(v1);
                    newtrivtx.push_back(v);
                    newtrivtx.push_back(v1);
                    newtrivtx.push_back(v2);
                    
                    vec3d nstart = (v1->pos()-b->pos()).cross(v2->pos()-b->pos());
                    vec3d nend   = (v1->pos()-newpos).cross(v2->pos()-newpos);
                    if( nstart.dot(nend) <= 0.0 )
                        return;
                    
                    if( edge_is_feature(b,v1) )
                        newfeat.push_back( halfedge(v,v1) );
                    if( edge_is_feature(b,v2) )
                        newfeat.push_back( halfedge(v,v2) );
                    if( edge_is_feature(a,v1) )
                        newfeat.push_back( halfedge(v,v1) );
                    if( edge_is_feature(a,v2) )
                        newfeat.push_back( halfedge(v,v2) );
                    
                    oldtrivtx.push_back(b);
                    oldtrivtx.push_back(v1);
                    oldtrivtx.push_back(v2);
                    remtri.push_back(t);
                }
            }
            for( int i=0; i<remtri.size(); i++ ){
                remove_triangle(remtri[i]);
            }
            bool failed=false;
            for( int i=0; i<newtrivtx.size(); i+=3 ){
                triangle* t = add_triangle( newtrivtx[i+0], newtrivtx[i+1], newtrivtx[i+2] );
                if( !t ){
                    failed=true;
                } else {
                    newtri.push_back( t );
                }
            }
            
            if( !failed ){
                if( (feat_edge && feat_edge_count < 6) || (!a->is_feature() && !b->is_feature()) )
                    a->pos() = (a->pos()+b->pos())/2.0;
                
                for( int i=0; i<newtrivtx.size(); i+=3 ){
                    for( int j=0; j<3; j++ ){
                        vertex* v0 = newtrivtx[i+j];
                        vertex* v1 = newtrivtx[i+(j+1)%3];
                        if( !tri_for_edge(v0,v1) || !tri_for_edge(v1,v0) )
                            fail_and_die("uh oh");
                    }
                }
                for( int i=0; i<newfeat.size(); i++ ){
                    mark_as_features( newfeat[i].first, newfeat[i].second );
                }
                remove_vertex(b);
            } else {
                //std::cout << "failed" << std::endl;
                for( int i=0; i<newtri.size(); i++ ){
                    remove_triangle(newtri[i]);
                }
                for( int i=0; i<oldtrivtx.size(); i+=3 ){
                    add_triangle( oldtrivtx[i+0], oldtrivtx[i+1], oldtrivtx[i+2] );
                }
            }
            return;
            
            /*
            std::vector<triangle*> dtris;
            std::vector<vertex*> oldtris;
            std::vector<vertex*> tris;
            std::vector<triangle*> newtris;
            for( vertex::triangle_iterator it=b->triangles_begin(); it!=b->triangles_end(); ++it ){
                triangle* t = *it;
                if( !t->has_vertex(a) ){
                    vertex *v1 = t->next(b);
                    vertex *v2 = t->next(v1);
                    tris.push_back(a);
                    tris.push_back(v1);
                    tris.push_back(v2);
                    dtris.push_back(t);
                    oldtris.push_back( t->vtx(0) );
                    oldtris.push_back( t->vtx(1) );
                    oldtris.push_back( t->vtx(2) );
                } else {
                    dtris.push_back(t);
                    oldtris.push_back( t->vtx(0) );
                    oldtris.push_back( t->vtx(1) );
                    oldtris.push_back( t->vtx(2) );
                }
            }

            // remove the old triangles
            for( int i=0; i<dtris.size(); i++ ){
                if( !check_mesh_pointers(dtris[i]) )
                    fail_and_die("mesh pointers incorrect");
                remove_triangle( dtris[i] );
            }
            
            // try to add the new triangles
            bool failed = false;
            for( int i=0; i<tris.size(); i+=3 ){
                vertex* v0 = tris[i+0];
                vertex* v1 = tris[i+1];
                vertex* v2 = tris[i+2];
                
                triangle* t = add_triangle(v0,v1,v2);
                if( t != NULL ){
                    newtris.push_back( t );
                } else {
                    failed = true;
                    break;
                }
            }
            //failed = true;
            if( failed ){
                for( int i=0; i<newtris.size(); i++ ){
                    remove_triangle( newtris[i] );
                }
                for( int i=0; i<oldtris.size(); i+=3 ){
                    add_triangle( oldtris[i+0], oldtris[i+1], oldtris[i+2] );
                }
                for( int i=0; i<oldtris.size(); i++ ){
                    if( !check_local_mesh( oldtris[i] ) )
                        fail_and_die("uhoh");
                }
                return;
            } else {
                for( int i=0; i<newtris.size(); i++ ){
                    triangle* t=newtris[i];
                    if( tri_for_edge( t->vtx(0),t->vtx(1) ) != t )
                        fail_and_die("uhoh");
                    if( tri_for_edge( t->vtx(1),t->vtx(2) ) != t )
                        fail_and_die("uhoh");
                    if( tri_for_edge( t->vtx(2),t->vtx(0) ) != t )
                        fail_and_die("uhoh");
                    if( !tri_for_edge( t->vtx(1),t->vtx(0) ) )
                        fail_and_die("uhoh");
                    if( !tri_for_edge( t->vtx(2),t->vtx(1) ) )
                        fail_and_die("uhoh");
                    if( !tri_for_edge( t->vtx(0),t->vtx(2) ) )
                        fail_and_die("uhoh");
                    if( !check_local_mesh( t->vtx(0) ) )
                        fail_and_die("uhoh");
                    if( !check_local_mesh( t->vtx(1) ) )
                        fail_and_die("uhoh");
                }
            }
            
            //remove_vertex(b);
            
            if( !check_local_mesh(a) )
                fail_and_die("error in local mesh");
            */
        }
        
        void split_edge( vertex* a, vertex* b ){
            triangle* T1 = tri_for_edge(a,b);
            triangle* T2 = tri_for_edge(b,a);
            if( !T1 || !T2 ) return;
            vertex *m = add_vertex( (a->pos()+b->pos())/2 );
            vertex *c = T1->other(a,b);
            vertex *d = T2->other(a,b);
            remove_triangle(T1);
            remove_triangle(T2);
            add_triangle( a, m, c );
            add_triangle( m, b, c );
            add_triangle( b, m, d );
            add_triangle( m, a, d );
            if( edge_is_feature( a, b ) ){
                mark_as_features( a, m );
                mark_as_features( m, b );
            }
        }
        
        void mark_as_features( vertex* a, vertex* b ){
            a->is_feature() = true;
            b->is_feature() = true;
            m_features.insert( halfedge(std::min(a,b),std::max(a,b)) );
        }
        
        bool edge_is_feature( vertex* a, vertex* b ){
            return m_features.find( halfedge(std::min(a,b),std::max(a,b)) ) != m_features.end();
        }

        void get_mesh( std::vector<double>& coord, std::vector<int>& tris ){
            int nextvid=0;
            std::map< vertex*, int > vid;
            for( mesh::vertex_iterator it=vertices_begin(); it!=vertices_end(); ++it ){
                vertex *v = *it;
                vid[v] = nextvid++;
                coord.push_back( v->pos()[0] );
                coord.push_back( v->pos()[1] );
                coord.push_back( v->pos()[2] );
            }
            for( mesh::triangle_iterator it=triangles_begin(); it!=triangles_end(); ++it ){
                triangle* t = *it;
                tris.push_back( vid[t->vtx(0)] );
                tris.push_back( vid[t->vtx(1)] );
                tris.push_back( vid[t->vtx(2)] );
            }
        }
        
        bool check_local_mesh( vertex* v ){
            bool ret=true;
            
            // first check that the vertex exists in the mesh
            if( !has_vertex(v) ){
                ret = false;
            }
            
            // now loop over the adjacent triangles
            for( vertex::triangle_iterator it=v->triangles_begin(); it!=v->triangles_end(); ++it ){
                triangle *t = *it;
                
                // check if the triangle exists in the mesh
                if( !has_triangle(t) ){
                    ret = false;
                }
                // check that the triangle actually uses the vertex
                if( !t->has_vertex(v) ){
                    ret = false;
                }
                
                // now loop over the edges
                for( int i=0; i<3; i++ ){
                    // check that the edge pointers point to the current triangle
                    if( tri_for_edge( t->vtx(i), t->vtx((i+1)%3)) != t ){
                        ret = false;
                    }
                    
                    // check that a neigboring triangle exists
                    if( !tri_for_edge( t->vtx((i+1)%3), t->vtx(i) ) ){
                        ret = false;
                    }
                }
            }
            return ret;
        }
        
        bool check_mesh(){
            for( vertex_iterator it=vertices_begin(); it!=vertices_end(); ++it ){
                if(!check_local_mesh(*it) ) return false;
            }
            return true;
        }
        
        bool check_mesh_pointers( vertex* v ){
            bool ret = true;
            for( vertex::triangle_iterator it=v->triangles_begin(); it!=v->triangles_end(); ++it ){
                triangle* t = *it;
                ret &= check_mesh_pointers( t );
            }
            return ret;
        }
        
        bool check_mesh_pointers( triangle* t ){
            if( tri_for_edge(t->vtx(0),t->vtx(1)) != t || tri_for_edge(t->vtx(1),t->vtx(2)) != t || tri_for_edge(t->vtx(2),t->vtx(0)) != t ){
                return false;
            }
            for( int i=0; i<3; i++ ){
                if( !t->vtx(0)->has_triangle(t) )
                    return false;
            }
            return true;
        }
        
        void fail_and_die( const char *msg ){
            std::vector<double> coord;
            std::vector<int> tri;
            get_mesh( coord, tri );
            save_obj( "remesher_dump.obj", coord, tri );
            std::cout << "FATAL ERROR: " << msg << std::endl;
            throw msg;
        }
        
        void fail_and_die( const char *msg, std::vector<double>& coord, std::vector<int>& tri ){
            save_obj( "remesher_error.obj", coord, tri );
            std::cout << "FATAL ERROR: " << msg << std::endl;
            throw msg;
        }
        
        void check_closed_manifold(){
            for( triangle_iterator it=triangles_begin(); it!=triangles_end(); ++it ){
                triangle* t = *it;
                if( !check_mesh_pointers(t) || !tri_for_edge(t->vtx(1), t->vtx(0)) || !tri_for_edge(t->vtx(2),t->vtx(1)) || !tri_for_edge(t->vtx(0),t->vtx(2)) ){
                    fail_and_die( "mesh is non_manifold" );
                }
            }
        }
        
    private:
        
        vertex_set     m_verts;
        triangle_set   m_tris;
        
        typedef std::pair<vertex*,vertex*> halfedge;
        std::map< halfedge, triangle*> m_tri_for_edge;
        
        std::set< halfedge > m_features;
    };
    
};

#endif