package picard;

import org.testng.Assert;
import org.testng.annotations.Test;


/**
 * Created by farjoun on 10/22/17.
 */
public class TestStringUtilSplitTest {

    @Test
    public void TestTestStringUtilSplit() {

        TestStringUtilSplit tester = new TestStringUtilSplit();
        Assert.assertTrue(tester.run() < 0.1, "We should stop using StringUtil.split as the StringTokenizer is 10% faster");
    }
}