#include <iostream>
using namespace std;

class Node {
public:
    double data;
    Node* next;

    Node(double val) {
        data = val;
        next = nullptr;
    }
};

class CircularLinkedList {
private:
    Node* head;

public:
    double invT; 
    int member_num; 
    CircularLinkedList(double val) {
        invT = val; 
        member_num = 0; 
        head = nullptr;
    }

    void insertAfter(Node* prevNode, double value) {
        if (prevNode == nullptr) return;

        Node* newNode = new Node(value);
        newNode->next = prevNode->next;
        prevNode->next = newNode;
    }

    void insertAtEnd(double value) {
        Node* newNode = new Node(value);
        if (head == nullptr) {
            head = newNode;
            head->next = head;
        } else {
            Node* temp = head;
            while (temp->next != head) {
                temp = temp->next;
            }
            temp->next = newNode;
            newNode->next = head;
        }
    }
    double findInterval(Node* prevNode){
        double val1  = prevNode->data;
        double val2  = prevNode->next->data;
        if (val2>val1){return val2-val1;}
        else{return invT-(val2-val1);} //Wrap around has occurred. 
    }

    void display() {
    if (head == nullptr) return;

    Node* temp = head;
    do {
        cout << temp->data << " ";
        temp = temp->next;
    } while (temp != head);
    cout << endl;
}
};


