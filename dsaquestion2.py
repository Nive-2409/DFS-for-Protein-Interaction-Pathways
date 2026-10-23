
# T9: DFS for Protein Interaction Pathways
# Modules:
# 1. Connected Components
# 2. Cycle Detection
# 3. Path Finding

class Stack:
    """Manually implemented stack for iterative DFS."""

    def __init__(self):
        self.items = []

    def push(self, item):
        self.items.append(item)

    def pop(self):
        if self.is_empty():
            return None
        return self.items.pop()

    def is_empty(self):
        return len(self.items) == 0


class ProteinNetwork:
    def __init__(self):
        self.proteins = []
        self.adjacency = {}

    def add_protein(self, protein):
        if protein not in self.adjacency:
            self.proteins.append(protein)
            self.adjacency[protein] = []

    def add_interaction(self, p1, p2):
        if p1 not in self.adjacency or p2 not in self.adjacency:
            raise ValueError("Both proteins must exist in the network.")

        if p1 == p2:
            raise ValueError("Self-interactions are not allowed.")

        if p2 not in self.adjacency[p1]:
            self.adjacency[p1].append(p2)
            self.adjacency[p2].append(p1)

    def display_network(self):
        print("\nProtein Interaction Network:")
        for protein in self.proteins:
            print(protein, "->", self.adjacency[protein])

    # MODULE 1: Connected Components
    def connected_components(self):
        visited = []
        components = []

        for start in self.proteins:
            if start not in visited:
                stack = Stack()
                stack.push(start)
                visited.append(start)
                component = []

                while not stack.is_empty():
                    current = stack.pop()
                    component.append(current)

                    for neighbour in self.adjacency[current]:
                        if neighbour not in visited:
                            visited.append(neighbour)
                            stack.push(neighbour)

                components.append(component)

        return components

    # MODULE 2: Cycle Detection
    def has_cycle(self):
        visited = []

        for start in self.proteins:
            if start in visited:
                continue

            stack = Stack()
            stack.push((start, None))
            visited.append(start)

            while not stack.is_empty():
                current, parent = stack.pop()

                for neighbour in self.adjacency[current]:
                    if neighbour not in visited:
                        visited.append(neighbour)
                        stack.push((neighbour, current))
                    elif neighbour != parent:
                        return True

        return False

    # MODULE 3: Path Finding
    def find_path(self, source, destination):
        if source not in self.adjacency:
            return None

        if destination not in self.adjacency:
            return None

        stack = Stack()
        visited = []
        parent = {}

        stack.push(source)
        visited.append(source)
        parent[source] = None

        while not stack.is_empty():
            current = stack.pop()

            if current == destination:
                path = []
                node = destination

                while node is not None:
                    path.append(node)
                    node = parent[node]

                path.reverse()
                return path

            for neighbour in self.adjacency[current]:
                if neighbour not in visited:
                    visited.append(neighbour)
                    parent[neighbour] = current
                    stack.push(neighbour)

        return []


def build_sample_network():
    network = ProteinNetwork()

    # 20 synthetic protein identifiers
    for i in range(1, 21):
        network.add_protein("P" + str(i))

    # Component 1: P1 to P8
    edges = [
        ("P1", "P2"),
        ("P1", "P3"),
        ("P2", "P4"),
        ("P3", "P4"),
        ("P4", "P5"),
        ("P5", "P6"),
        ("P6", "P7"),
        ("P7", "P8"),

        # Component 2: P9 to P14
        ("P9", "P10"),
        ("P10", "P11"),
        ("P11", "P12"),
        ("P12", "P13"),
        ("P13", "P14"),

        # Component 3: P15 to P18
        ("P15", "P16"),
        ("P16", "P17"),
        ("P17", "P18"),
        ("P18", "P15"),

        # Component 4: P19 and P20
        ("P19", "P20")
    ]

    for p1, p2 in edges:
        network.add_interaction(p1, p2)

    return network


def main():
    network = build_sample_network()

    while True:
        print("\n===== PROTEIN INTERACTION PATHWAYS =====")
        print("1. Display Network")
        print("2. Find Connected Components")
        print("3. Detect Cycle")
        print("4. Find Path Between Proteins")
        print("5. Exit")

        choice = input("Enter your choice: ").strip()

        if choice == "1":
            network.display_network()

        elif choice == "2":
            components = network.connected_components()
            print("\nConnected Components:")
            for i, component in enumerate(components, start=1):
                print("Component", i, ":", component)
            print("Total components:", len(components))

        elif choice == "3":
            if network.has_cycle():
                print("Cycle detected in the protein network.")
            else:
                print("No cycle detected.")

        elif choice == "4":
            source = input("Enter source protein (e.g., P1): ").strip().upper()
            destination = input(
                "Enter destination protein (e.g., P6): "
            ).strip().upper()

            if source not in network.adjacency:
                print("Error: Source protein does not exist.")
            elif destination not in network.adjacency:
                print("Error: Destination protein does not exist.")
            else:
                path = network.find_path(source, destination)

                if path:
                    print("Path:", " -> ".join(path))
                else:
                    print("No path exists between the proteins.")

        elif choice == "5":
            print("Exiting the program.")
            break

        else:
            print("Invalid choice. Enter a number from 1 to 5.")


if __name__ == "__main__":
    main()
