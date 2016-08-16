// build this guy with:
//  g++ -o test_graph.exe main.cc Graph.cc -I/users/plato/rosetta/SVNROSETTA/mini_branch/mini/src
// ../../utility/pointer/ReferenceCount.cc
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/graph/Graph.hh>

#include <iostream>

using namespace core::graph;

void print_edges( Node const * node );

int main()
{

	std::cout << "test" << std::endl;

	GraphOP g = new Graph( 10 );
	g->add_edge(1,2);
	g->add_edge(2,3);
	g->add_edge(3,4);
	g->add_edge(4,5);
	g->add_edge(5,6);
	g->add_edge(6,7);
	g->add_edge(6,10);
	g->add_edge(4,6);
	g->add_edge(7,8);
	g->add_edge(8,9);
	g->add_edge(9,10);

	g->print_vertices();

	std::cout << "size of g: " << g->getTotalMemoryUsage() << " bytes." << std::endl;

	Node * node6 = g->get_node( 6 );

	std::cout << "edges for node 6" << std::endl;

	print_edges( node6 );

	Edge* edge_6_7 = g->find_edge(6,7);
	g->delete_edge( edge_6_7 ); edge_6_7 = 0;

	print_edges( node6 );

	Edge* edge_6_10 = g->find_edge(6,10);
	g->delete_edge( edge_6_10 );

	print_edges( node6 );

	g->print_vertices();
	std::cout << "size of g: " << g->getTotalMemoryUsage() << " bytes." << std::endl;

	//as return is called, g passes out of scope and the Graph is deallocated

	GraphOP g2 = new Graph();
	g2->set_num_nodes( 10 );
	g2->add_edge(1,3);
	g2->add_edge(2,4);
	g2->add_edge(3,5);
	g2->add_edge(4,6);
	g2->add_edge(5,7);
	g2->add_edge(6,8);
	g2->add_edge(6,10);
	g2->add_edge(4,7);
	g2->add_edge(7,8);
	g2->add_edge(8,10);
	g2->add_edge(9,10);

	std::cout << "g2" << std::endl;
	g2->print_vertices();

	std::cout << "*g2 = *g" << std::endl;
	*g2 = *g;

	std::cout << "g2" << std::endl;
	g2->print_vertices();


	std::cout << "Adding Loops  to g2" << std::endl;
	g2->add_edge(1,1);
	g2->add_edge(7,7);
	print_edges( g2->get_node( 7 ) );
	Edge* edge_7_7 = g2->find_edge( 7,7 );
	g->delete_edge( edge_7_7 );
				print_edges( g2->get_node( 7 ) );
	g2->print_vertices();
	return 1;

}

void print_edges( Node const * node )
{
	std::cout << "All edges for node: " << node->get_node_index() << std::endl;
	for( Node::EdgeListConstIter
		citer = node->const_edge_list_begin();
		citer != node->const_edge_list_end(); ++citer )
	{
		std::cout << "Edge between " << (*citer)->get_first_node_ind() << " ";
		std::cout << (*citer)->get_second_node_ind() << std::endl;
	}

	std::cout << "Lower edges for node: " << node->get_node_index() << std::endl;
	for( Node::EdgeListConstIter
		citer = node->const_lower_edge_list_begin();
		citer != node->const_lower_edge_list_end(); ++citer )
	{
		std::cout << "Edge between " << (*citer)->get_first_node_ind() << " ";
		std::cout << (*citer)->get_second_node_ind() << std::endl;
	}

	std::cout << "Upper edges for node: " << node->get_node_index() << std::endl;
	for( Node::EdgeListConstIter
		citer = node->const_upper_edge_list_begin();
		citer != node->const_upper_edge_list_end(); ++citer )
	{
		std::cout << "Edge between " << (*citer)->get_first_node_ind() << " ";
		std::cout << (*citer)->get_second_node_ind() << std::endl;
	}

}
