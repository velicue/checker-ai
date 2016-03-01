#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <vector>
#include "rlutil.h"
#include <algorithm>

using namespace std;

const int MAXN = 9;

class Point
{
public:
    int x, y;
    Point() {}
    Point(const Point& p) {x = p.x; y = p.y;}
    Point(int x_, int y_) {x = x_; y = y_;}
    bool operator == (const Point &a) const {return x == a.x && y == a.y;}
};

class resultNode
{
public:
    vector<Point> nextmove;
    int prescore;
    int evaluation;
    bool operator < (const resultNode &a) const {return prescore > a.prescore;}
};

int dx[6] = {-1, -1, 0, 0, 1, 1};
int dy[6] = {0, 1, -1, 1, -1, -0};

int board[MAXN][MAXN];
const int initscore = 2 * (MAXN - 1) * 1 + 2 * (MAXN - 2) * 2 + 2 * (MAXN - 3) * 3 + 2 * (MAXN - 4) * 4;
const int initsidescore = (MAXN - 1) * 4 + (MAXN - 3) * 3 + (MAXN - 5) * 2 + (MAXN - 7) * 1;
int score[8] = {initscore, initscore, initsidescore, initsidescore, initsidescore, initsidescore, 10, 10};

void printboard(int board[MAXN][MAXN], vector<Point> validlist = vector<Point>(), vector<Point> movelist = vector<Point>())
{
    int ct[2] = {0, 0};

    for (int i = 0; i < MAXN * 2 - 1; i++)
    {
        for (int j = 0; j < abs(MAXN - i - 1); j++)
            printf(" ");
        for (int j = MAXN - 1; j >= 0; j--)
        {
            if (j > i) continue;
            if (i - j >= MAXN) continue;

            printf(" ");

            int foundmove = -1;
            for (int k = 0; k < movelist.size(); k++)
            {
                if (j == movelist[k].x && i - j == movelist[k].y)
                {
                    foundmove = k;
                    break;
                }
            }
            if (foundmove != -1)
            {
                if (board[j][i - j] == 0) rlutil::setColor(rlutil::YELLOW);
                else if (board[j][i - j] == 1) rlutil::setColor(rlutil::LIGHTGREEN);
                else rlutil::setColor(rlutil::LIGHTRED);

                printf("%c", 'a' + foundmove);
                rlutil::setColor(rlutil::GREY);
                continue;
            }

            if (board[j][i - j] == -1)
            {
                int index = -1;
                for (int k = 0; k < validlist.size(); k++)
                {
                    if (validlist[k] == Point(j, i - j))
                    {
                        index = k;
                        break;
                    }
                }
                if (index == -1)
                    printf(".");
                else
                    printf("%c", 'a' + index);
                continue;
            }
            if (board[j][i - j] == 0)
                rlutil::setColor(rlutil::YELLOW);
            if (board[j][i - j] == 1)
                rlutil::setColor(rlutil::LIGHTGREEN);
            printf("%d", ct[board[j][i - j]]++);
            rlutil::setColor(rlutil::GREY);
        }
        printf("\n");
    }
}

bool check(int board[MAXN][MAXN], Point p, int player)
{
    if (p.x < 0 || p.x >= MAXN || p.y < 0 || p.y >= MAXN || (board[p.x][p.y] != player)) return false;
    return true;
}

void domove(int board[MAXN][MAXN], int score[2], vector<Point> plist, int player, bool rev = false)
{
    int x0, y0, xn, yn;
    if (!rev)
    {
        x0 = plist[0].x;
        y0 = plist[0].y;
        xn = plist[plist.size() - 1].x;
        yn = plist[plist.size() - 1].y;
    } else {
        xn = plist[0].x;
        yn = plist[0].y;
        x0 = plist[plist.size() - 1].x;
        y0 = plist[plist.size() - 1].y;
    }

    board[x0][y0] = -1;
    board[xn][yn] = player;
    if (player == 0)
    {
        score[player] += xn + yn - x0 - y0;
        score[player + 2] += xn - x0;
        score[player + 4] += yn - y0;

        if (x0 + y0 >= MAXN * 2 - 5 && xn + yn < MAXN * 2 - 5)
            score[player + 6] -= 1;
        if (xn + yn >= MAXN * 2 - 5 && x0 + y0 < MAXN * 2 - 5)
            score[player + 6] += 1;
    } else {
        score[player] -= xn + yn - x0 - y0;
        score[player + 2] -= xn - x0;
        score[player + 4] -= yn - y0;
        if (x0 + y0 < 4 && xn + yn >= 4)
            score[player + 6] -= 1;
        if (xn + yn < 4 && x0 + y0 >= 4)
            score[player + 6] += 1;
    }
}

Point findmove(int board[MAXN][MAXN], int no, int player)
{
    int ct = 0;
    for (int i = 0; i < MAXN * 2 - 1; i++)
    {
        for (int j = MAXN - 1; j >= 0; j--)
        {
            if (j > i) continue;
            if (i - j >= MAXN) continue;
            if (board[j][i - j] == player)
            {
                if (no == ct++)
                {
                    return Point(j, i - j);
                }
            }
        }
    }
    return Point(-1, -1);
}


const double coef = 100000.0;
const int INF = 2147483647;

vector<Point> validmoves(int board[MAXN][MAXN], Point tp, int player, int mode)
{
    vector<Point> result;

    if (mode == 0)
    {
        for (int i = 0; i < 6; i++)
        {
            Point np = Point(tp.x + dx[i], tp.y + dy[i]);
            if (!check(board, np, -1)) continue;
            result.push_back(np);
        }
    }
    for (int i = 0; i < 6; i++)
    {
        Point np = Point(tp.x + dx[i], tp.y + dy[i]);
        int k = 1;
        while (check(board, np, -1))
        {
            np.x += dx[i];
            np.y += dy[i];
            k++;
        }
        if (!check(board, np, 0) && !check(board, np, 1)) continue;
        while (k--)
        {
            np.x += dx[i];
            np.y += dy[i];
            if (!check(board, np, -1)) break;
        }
        if (!check(board, np, -1)) continue;
        result.push_back(np);
    }
    return result;
}

vector<Point> temp;
vector<vector<Point>> tempmoves;

boolean visited[MAXN][MAXN];

void generatemoves(int board[MAXN][MAXN], Point tp, int player, int mode)
{
    if (visited[tp.x][tp.y]) return;
    visited[tp.x][tp.y] = true;
    temp.push_back(tp);
    if (temp.size() > 1) tempmoves.push_back(temp);
    if (mode == -1)
    {
        temp.pop_back();
        return;
    }
    vector<Point> validlist = validmoves(board, tp, player, mode);
    for (int i = 0; i < validlist.size(); i++)
    {
        bool found = false;
        for (int j = 0; j < 6; j++)
        {
            if (tp.x + dx[j] == validlist[i].x && tp.y + dy[j] == validlist[i].y)
            {
                found = true;
                break;
            }
        }
        if (found)
            generatemoves(board, validlist[i], player, -1);
        else
            generatemoves(board, validlist[i], player, 1);
    }
    temp.pop_back();
}

int tempboard[MAXN][MAXN];
int tempscore[8];

int R = 2;

class HashNode
{
public:
    int board[MAXN][MAXN];
    int player;
    resultNode res;
    bool used;
    int dep;
    int history;
};

const int MAXH = 4999997;
HashNode hashtable[MAXH];
int HISTORY = 0;

int calchash(int board[MAXN][MAXN], int player)
{
    int sum = player;
    for (int i = 0; i < MAXN; i++)
    {
        for (int j = 0; j < MAXN; j++)
        {
            sum *= 133;
            sum += board[i][j];
        }
    }
    while (sum < 0) sum += MAXH;
    sum %= MAXH;
    return sum;
}

int gethash(int board[MAXN][MAXN], int player, int dep)
{
    int hc = calchash(board, player);
    if (!hashtable[hc].used) return MAXH;
    if (hashtable[hc].player != player) return MAXH;
    if (hashtable[hc].dep + 1 < dep) return MAXH;
    for (int i = 0; i < MAXN; i++)
    {
        for (int j = 0; j < MAXN; j++)
        {
            if (board[i][j] != hashtable[hc].board[i][j])
            {
                return MAXH;
            }
        }
    }

    return hc;
}

void storehash(int board[MAXN][MAXN], int player, int dep, resultNode res)
{
    int hc = calchash(board, player);
    if (hashtable[hc].used = true && hashtable[hc].dep >= dep && hashtable[hc].history >= HISTORY - 1)
    {
        return;
    }
    hashtable[hc].used = true;
    hashtable[hc].dep = dep;
    hashtable[hc].player = player;
    hashtable[hc].res = res;
    memcpy(hashtable[hc].board, board, sizeof(int) * MAXN * MAXN);
    hashtable[hc].history = HISTORY;
}

int STARTTIME = 0;
int TIMELIMIT = 60;
int MAXDEPTH = 12;
int MAXBRANCH = 18;

resultNode DFS(int dep, int player, int alpha, int beta, int branch, int allownull = 0)
{
    /*int hc = gethash(tempboard, player, dep);
    if (hc != MAXH && hc >= 0)
    {
        printf("HIT!\n");
        return hashtable[hc].res;
    }
    */

    resultNode res;
    res.evaluation = alpha;

    //printf("DEP: %d\n", dep);


    if (tempscore[player] == 0) return res.evaluation = INF - 2 - MAXDEPTH + dep, res;
    if (tempscore[1 - player] == 0) return res.evaluation = -(INF - 2 - MAXDEPTH + dep), res;
    if (dep == 0)
    {
        double t1 = max(max(tempscore[1 - player], tempscore[1 - player + 2]), tempscore[1 - player + 4]);
        t1 += (tempscore[1 - player + 6] * (120 / t1) / 2.0);
        double t2 = max(max(tempscore[player], tempscore[player + 2]), tempscore[player + 4]);
        t2 += (tempscore[player + 6] * (120 / t2) / 2.0);
        return res.evaluation = int(coef * (log(t1) - log(t2 - t2 / 60.0))), res;
    }

    //NULL MOVE
    if (dep > R + 1 && allownull < 2)
    {
        int val = -DFS(dep - 1 - R, 1 - player, -beta, -beta + 1, branch * 4 / 5, allownull + 1).evaluation;
        if (val >= beta)
        {
            return res.evaluation = val, res;
        }
    }

    vector<resultNode> currentmoves;

    tempmoves.clear();
    for (int i = 0; i < MAXN; i++)
    {
        for (int j = 0; j < MAXN; j++)
        {
            if (tempboard[i][j] != player) continue;
            temp.clear();
            memset(visited, 0, sizeof(visited));
            tempboard[i][j] = -1;
            generatemoves(tempboard, Point(i, j), player, 0);
            tempboard[i][j] = player;
        }
    }


    for (int i = 0; i < tempmoves.size(); i++)
    {
        resultNode currentmove;
        currentmove.nextmove = tempmoves[i];
        int len = tempmoves[i].size();
        int currentscore = max(max(tempscore[player], tempscore[player + 2]), tempscore[player + 4]);
        int score0 = tempmoves[i][len - 1].x + tempmoves[i][len - 1].y - tempmoves[i][0].x - tempmoves[i][0].y;
        int score1 = tempmoves[i][len - 1].x - tempmoves[i][0].x;
        int score2 = tempmoves[i][len - 1].y - tempmoves[i][0].y;

        if (player == 0)
        {
            score0 = -score0;
            score1 = -score1;
            score2 = -score2;
        }

        int newscore = max(max(tempscore[player] - score0, tempscore[player + 2] - score1), tempscore[player + 4] - score2);
        currentmove.prescore = currentscore - newscore;

        currentmoves.push_back(currentmove);
    }
    sort(currentmoves.begin(), currentmoves.end());

    int b = beta;
    int foundPV = false;


    for (int i = 0; i < currentmoves.size(); i++)
    {
        domove(tempboard, tempscore, currentmoves[i].nextmove, player);
        int nextdepth = dep - 1;
        if (foundPV && nextdepth > 0) nextdepth -= max(nextdepth / 3, 1);
        if (i > branch && nextdepth > 0) nextdepth -= max(nextdepth / 4, 1);
        if (i > 2 * branch && nextdepth > 0) nextdepth -= max(nextdepth / 4, 1);
        if (i > 3 * branch && nextdepth > 0) nextdepth -= (nextdepth + 2) / 4;


        int val = -DFS(nextdepth, 1 - player, -b, -res.evaluation, branch * 4 / 5).evaluation;
        if (val > res.evaluation && val < beta)
        {
            if (i > 0) val = -DFS(dep - 1, 1 - player, -beta, -res.evaluation, branch * 4 / 5).evaluation;
            if (i > 0 || dep != MAXDEPTH) foundPV = true;
        }
        if (val > res.evaluation)
        {
            res = currentmoves[i];
            res.evaluation = val;
        }
        domove(tempboard, tempscore, currentmoves[i].nextmove, player, true);
        if (res.evaluation >= beta) return res;
        b = res.evaluation + 1;
        if (time(NULL) - STARTTIME > TIMELIMIT) return res;
    }
    if (dep == MAXDEPTH && res.nextmove.size() == 0)
    {
        printf("BUG\n");
        printf("%d\n", currentmoves.size());
    }

    //storehash(board, player, dep, res);
    return res;
}


vector<Point> cgen(int player)
{
    int t1 = time(NULL);
    STARTTIME = t1;
    memcpy(tempboard, board, sizeof(board));
    memcpy(tempscore, score, sizeof(score));
    resultNode res = DFS(MAXDEPTH, player, -INF, INF, MAXBRANCH);
    printf("TIME COSTED: %d\n", time(NULL) - t1);
    printf("CONFIDENCE: %d\n", res.evaluation);
    return res.nextmove;
}

int main()
{
    srand(time(NULL));
    memset(board, -1, sizeof(board));
    int currentPlayer = 0;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4 - i; j++)
        {
            board[i][j] = 1;
        }
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4 - i; j++)
        {
            board[MAXN - i - 1][MAXN - j - 1] = 0;
        }
    }

    printboard(board);

    int win = -1;
    int step = 0;
    HISTORY = 0;
    while (win == -1)
    {
        HISTORY++;
        if (currentPlayer == 0) rlutil::setColor(rlutil::YELLOW); else rlutil::setColor(rlutil::LIGHTGREEN);
        if (currentPlayer < 2) printf("C: \n"); else printf("U: \n");
        rlutil::setColor(rlutil::GREY);
        vector<Point> plist;
        if (currentPlayer < 2)
        {
            plist = cgen(currentPlayer);
            for (int i = 0; i < plist.size(); i++)
            {
                printf("%d %d\n", plist[i].x, plist[i].y);
            }
        } else {

            vector<Point> validlist;
            int no;
            scanf("%d", &no);
            Point p;
            while (p = findmove(board, no, currentPlayer), validlist = validmoves(board, p, currentPlayer, 0), validlist.size() == 0)
            {
                printf("Reinput: \n");
                scanf("%d", &no);
            }
            plist.push_back(p);

            /*
            printf("TARGET:\n");
            int tx = -1, ty = -1;
            while (!check(board, Point(tx, ty), -1))
                scanf("%d%d", &tx, &ty);
            board[p.x][p.y] = -1;
            board[tx][ty] = currentPlayer;
            plist.push_back(Point(tx, ty));
            */
            board[p.x][p.y] = -1;
            while (true)
            {
                printboard(board, validlist);
                char c = 0;
                scanf("%c", &c);
                while (c == '\n' || c - 'a' >= (int)validlist.size() || ((c < 'a' || c > 'z') && plist.size() == 1))
                    scanf("%c", &c);
                if (c < 'a' || c > 'z') break;
                int t = c - 'a';
                plist.push_back(validlist[t]);
                if (plist.size() == 2)
                {
                    bool found = false;
                    for (int i = 0; i < 6; i++)
                    {
                        if (plist[1].x == plist[0].x + dx[i] && plist[1].y == plist[0].y + dy[i])
                        {
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                }
                validlist = validmoves(board, validlist[t], currentPlayer, 1);
                if (validlist.size() == 0) break;
            }
            board[p.x][p.y] = currentPlayer;

        }
        domove(board, score, plist, currentPlayer);
        printboard(board, vector<Point>(), plist);
        printf("%d %d\n", max(max(score[0], score[2]), score[4]), max(max(score[1], score[3]), score[5]));
        if (score[0] == 0) win = 0;
        if (score[1] == 0) win = 1;
        currentPlayer = 1 - currentPlayer;
        step++;
    }
    if (win == 0) printf("CWINS!\n");
    else printf("UWINS!\n");

    return 0;
}
